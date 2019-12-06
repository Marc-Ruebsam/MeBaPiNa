checkpoint guppy:
    input:
        "00_rawdata/{run}/fast5"
    output:
        list(filter(None,[
        directory("01_processeddata/{run}/basecall/pass"),
        ## evaluation of Lambda calibration strands only when specified
        (config["guppy"]["lam_DCS"] and
        directory("01_processeddata/{run}/basecall/calibration_strands")),
        "01_processeddata/{run}/basecall/sequencing_summary.txt"
        ]))
    log:
        "01_processeddata/{run}/basecall/MeBaPiNa_guppy.log"
    benchmark:
        "01_processeddata/{run}/basecall/MeBaPiNa_guppy.benchmark.tsv"
    version:
        subprocess.check_output("guppy_basecaller --version | awk '{print $NF}'", shell=True)
    threads:
        8
    params:
        ("--flowcell " + config["guppy"]["flowcell"]),
        ("--kit " + config["guppy"]["seq_kit"]),
        ("--calib_detect" if config["guppy"]["lam_DCS"] else ""), ## includes detection of lambda clibration strands
        "--barcode_kits " + BAC_KIT, ## always includes demultiplexing (all reads marked as unclassified if no barcodes were used)
        # "--require_barcodes_both_ends",
        "--qscore_filtering",
        ("--min_qscore " + config["guppy"]["q_cut"]),
        ("--device cuda:all:100% "
        + "--gpu_runners_per_device 6 " 
        + "--chunks_per_runner 1536 " 
        + "--chunk_size 1000 " 
        + "--chunks_per_caller 10000 " 
        + "--num_barcode_threads 8" 
        if config["machine"]["gpu"] else ""),
        "--compress_fastq",
        "--fast5_out"
    shell:
        "guppy_basecaller --num_callers {threads} {params} "
        "--save_path 01_processeddata/{wildcards.run}/basecall "
        "--input_path {input} > {log} 2>&1"
