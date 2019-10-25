rule guppy:
    input:
        "00_rawdata/{run}/fast5"
    output:
        # directory("01_processeddata/{run}/workspace"),
        # directory("01_processeddata/{run}/pass"),
        # directory("01_processeddata/{run}/calibration_strands"),
        "01_processeddata/{run}/basecall/sequencing_summary.txt" ## used as dummy for the other folders
    log:
        "01_processeddata/{run}/basecall/MeBaPiNa_guppy.log"
    benchmark:
        "01_processeddata/{run}/basecall/MeBaPiNa_guppy.benchmark.tsv"
    version:
        subprocess.check_output("guppy_basecaller --version | awk '{print $NF}'", shell=True)
    threads:
        config["machine"]["cpu"]
    params:
        ("--flowcell " + config["guppy"]["flowcell"]),
        ("--kit " + config["guppy"]["seq_kit"]),
        ("--calib_detect" if config["guppy"]["lam_DCS"] else ""), ## includes detection of lambda clibration strands
        ("--barcode_kits " + config["guppy"]["bac_kit"] if config["guppy"]["bac_kit"] else ""), ## includes demultiplexing
        "--qscore_filtering",
        ("--min_qscore " + config["guppy"]["q_cut"]),
        ("--device cuda:all:100%" if config["machine"]["gpu"] else ""),
        "--compress_fastq",
        "--fast5_out"
    shell:
        "guppy_basecaller --num_callers {threads} {params} "
        "--save_path 01_processeddata/{run}/basecall "
        "--input_path {input} > {log} 2>&1"
