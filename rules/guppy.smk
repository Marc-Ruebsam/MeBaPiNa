rule guppy:
    input:
        "00_rawdata/{run}/fast5"
    output:
        temp(directory("01_basecalling/{run}/workspace")),
        protected(directory("01_basecalling/{run}/pass")),
        protected(directory("01_basecalling/{run}/fail")),
        protected(directory("01_basecalling/{run}/calibration_strands")),
        protected("01_basecalling/{run}/sequencing_summary.txt")
    log:
        "01_basecalling/{run}/MeBaPiNa.log"
    benchmark:
        "01_basecalling/{run}/MeBaPiNa.benchmark.log"
    threads:
        config["machine"]["cpu"]
    params: 
        ("--flowcell " + config["guppy"]["flowcell"]),
        ("--kit " + config["guppy"]["seq_kit"]),
        (config["guppy"]["lam_DCS"] and 
        "--calib_detect"),
        (config["guppy"]["bac_kit"] and 
        "--barcode_kits " + config["guppy"]["bac_kit"]),
        ("--qscore_filtering"),
        ("--min_qscore " + config["guppy"]["q_cut"]),
        (config["machine"]["gpu"] and 
        "--device cuda:all:100%"),
        (lambda wildcards, threads: "--num_callers " + str(threads)), ## have to use function; function always has to have wildcards argument first
        ("--compress_fastq"),
        ("--fast5_out"),
        ("--save_path 01_basecalling/{run}")
    shell:
        "guppy_basecaller {params} --input_path {input} / 2> {log}"
