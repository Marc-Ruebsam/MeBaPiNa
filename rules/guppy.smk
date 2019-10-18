rule guppy:
    input:
        "00_rawdata/{run}/fast5"
    output:
        directory("01_basecalling/{run}/workspace"),
        directory("01_basecalling/{run}/pass"),
        "01_basecalling/{run}/sequencing_summary.txt"
    log:
        "01_basecalling/{run}/MeBaPiNa.log"
    benchmark:
        "01_basecalling/{run}/MeBaPiNa.benchmark.tsv"
    version:
        subprocess.check_output("guppy_basecaller --version | awk '{print $NF}'", shell=True)
    threads:
        config["machine"]["cpu"]
    params: 
        ("--flowcell " + config["guppy"]["flowcell"]),
        ("--kit " + config["guppy"]["seq_kit"]),
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

rule guppy_lam:
    input:
        "00_rawdata/{run}/fast5"
    output:
        directory("01_basecalling/{run}/workspace"),
        directory("01_basecalling/{run}/pass"),
        directory("01_basecalling/{run}/calibration_strands"),
        "01_basecalling/{run}/sequencing_summary.txt"
    log:
        "01_basecalling/{run}/MeBaPiNa.log"
    benchmark:
        "01_basecalling/{run}/MeBaPiNa.benchmark.tsv"
    version:
        subprocess.check_output("guppy_basecaller --version | awk '{print $NF}'", shell=True)
    threads:
        config["machine"]["cpu"]
    params: 
        ("--flowcell " + config["guppy"]["flowcell"]),
        ("--kit " + config["guppy"]["seq_kit"]),
        ("--calib_detect"), ## includes detection of lambda clibration strands
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

rule guppy_bac:
    input:
        "00_rawdata/{run}/fast5"
    output:
        directory("01_basecalling/{run}/workspace"),
        directory("01_basecalling/{run}/pass"),
        "01_basecalling/{run}/sequencing_summary.txt"
    log:
        "01_basecalling/{run}/MeBaPiNa.log"
    benchmark:
        "01_basecalling/{run}/MeBaPiNa.benchmark.tsv"
    version:
        subprocess.check_output("guppy_basecaller --version | awk '{print $NF}'", shell=True)
    threads:
        config["machine"]["cpu"]
    params: 
        ("--flowcell " + config["guppy"]["flowcell"]),
        ("--kit " + config["guppy"]["seq_kit"]),
        ("--barcode_kits " + config["guppy"]["bac_kit"]), ## includes demultiplexing
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

rule guppy_lam_bac:
    input:
        "00_rawdata/{run}/fast5"
    output:
        directory("01_basecalling/{run}/workspace"),
        directory("01_basecalling/{run}/pass"),
        directory("01_basecalling/{run}/calibration_strands"),
        "01_basecalling/{run}/sequencing_summary.txt"
    log:
        "01_basecalling/{run}/MeBaPiNa.log"
    benchmark:
        "01_basecalling/{run}/MeBaPiNa.benchmark.tsv"
    version:
        subprocess.check_output("guppy_basecaller --version | awk '{print $NF}'", shell=True)
    threads:
        config["machine"]["cpu"]
    params: 
        ("--flowcell " + config["guppy"]["flowcell"]),
        ("--kit " + config["guppy"]["seq_kit"]),
        ("--calib_detect"), ## includes detection of lambda clibration strands
        ("--barcode_kits " + config["guppy"]["bac_kit"]), ## includes demultiplexing
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

ruleorder: guppy > guppy_lam > guppy_bac > guppy_lam_bac
