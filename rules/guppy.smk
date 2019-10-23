rule guppy:
    input:
        "00_rawdata/{run}/fast5"
    output:
        directory("01_processeddata/{run}/workspace"),
        "01_processeddata/{run}/pass/{reads}.fastq.gz",
        "01_processeddata/{run}/sequencing_summary.txt"
    log:
        "01_processeddata/{run}/MeBaPiNa_guppy.log"
    benchmark:
        "01_processeddata/{run}/MeBaPiNa_guppy.benchmark.tsv"
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
        ("--compress_fastq"),
        ("--fast5_out"),
        ("--save_path 01_processeddata/{run}")
    shell:
        "guppy_basecaller --num_callers {threads} {params} --input_path {input} / 2> {log}"

rule guppy_lam:
    input:
        "00_rawdata/{run}/fast5"
    output:
        directory("01_processeddata/{run}/workspace"),
        "01_processeddata/{run}/pass/{reads}.fastq.gz",
        directory("01_processeddata/{run}/calibration_strands"),
        "01_processeddata/{run}/sequencing_summary.txt"
    log:
        "01_processeddata/{run}/MeBaPiNa_guppy_lam.log"
    benchmark:
        "01_processeddata/{run}/MeBaPiNa_guppy_lam.benchmark.tsv"
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
        ("--compress_fastq"),
        ("--fast5_out"),
        ("--save_path 01_processeddata/{run}")
    shell:
        "guppy_basecaller --num_callers {threads} {params} --input_path {input} / 2> {log}"

rule guppy_bac:
    input:
        "00_rawdata/{run}/fast5"
    output:
        directory("01_processeddata/{run}/workspace"),
        "01_processeddata/{run}/pass/{reads}.fastq.gz",
        "01_processeddata/{run}/sequencing_summary.txt"
    log:
        "01_processeddata/{run}/MeBaPiNa_guppy_bac.log"
    benchmark:
        "01_processeddata/{run}/MeBaPiNa_guppy_bac.benchmark.tsv"
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
        ("--compress_fastq"),
        ("--fast5_out"),
        ("--save_path 01_processeddata/{run}")
    shell:
        "guppy_basecaller --num_callers {threads} {params} --input_path {input} / 2> {log}"

rule guppy_lam_bac:
    input:
        "00_rawdata/{run}/fast5"
    output:
        directory("01_processeddata/{run}/workspace"),
        "01_processeddata/{run}/pass/{reads}.fastq.gz",
        directory("01_processeddata/{run}/calibration_strands"),
        "01_processeddata/{run}/sequencing_summary.txt"
    log:
        "01_processeddata/{run}/MeBaPiNa_guppy_lam_bac.log"
    benchmark:
        "01_processeddata/{run}/MeBaPiNa_guppy_lam_bac.benchmark.tsv"
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
        ("--compress_fastq"),
        ("--fast5_out"),
        ("--save_path 01_processeddata/{run}")
    shell:
        "guppy_basecaller --num_callers {threads} {params} --input_path {input} / 2> {log}"

ruleorder: guppy > guppy_lam > guppy_bac > guppy_lam_bac
