#################
## BASECALLING ##
#################

## BASECALL DEMULTIPLEX ##
##########################

rule move_raw:
    output:
        directory("{tmp}00_raw_data/{run}/fast5")
    shell:
        "long_path=$(find {wildcards.tmp}00_raw_data/ -name \"{wildcards.run}\"); "
        "mv ${{long_path}} \"{wildcards.tmp}00_raw_data/{wildcards.run}\"; "
        "find {wildcards.tmp}00_raw_data/ -depth -type d -empty -delete" #!# delete empty directories

checkpoint basecall_raw:
    input:
        "{tmp}00_raw_data/{run}/fast5"
    output:
        ([directory("{tmp}01_processed_data/01_basecalling/{run}/pass"),
        "{tmp}01_processed_data/01_basecalling/{run}/sequencing_summary.txt"] +
        [x for x in [directory("{tmp}01_processed_data/01_basecalling/{run}/calibration_strands")] if LAM_DCS ]) ## evaluation of Lambda calibration strands only when specified
    log:
        "{tmp}01_processed_data/01_basecalling/{run}/MeBaPiNa_basecall_raw.log"
    benchmark:
        "{tmp}01_processed_data/01_basecalling/{run}/MeBaPiNa_basecall_raw.benchmark.tsv"
    version:
        subprocess.check_output("guppy_basecaller --version | awk '{print $NF}'", shell=True)
    threads:
        22
    params:
        "--flowcell " + FLOWCELL,
        "--kit " + SEQ_KIT,
        "--barcode_kits " + BAC_KIT, ## always includes demultiplexing (all reads marked as unclassified if no barcodes were used)
        ("--calib_detect" if LAM_DCS else ""), ## includes detection of lambda clibration strands
        "--min_score 70", ## minimum score for barcode match (default 60)
        # "--require_barcodes_both_ends", ## very stringent
        # "--min_score_rear_override 70"
        # "--detect_mid_strand_barcodes", ## very slow
        # "--min_score_mid_barcodes 70"
        # "--trim_barcodes",
        # "--num_extra_bases_trim 10",
        "--qscore_filtering",
        ("--min_qscore " + config["filtering"]["q_min"]),
        ("--device cuda:all:100% "
        + "--gpu_runners_per_device 6 "
        + "--chunks_per_runner 1536 "
        + "--chunk_size 1000 "
        + "--chunks_per_caller 10000 "
        + "--num_barcode_threads 8"
        if config["workstation"]["gpu"] else ""),
        # "--compress_fastq", #!# if commented does not compress fastq, because of downstream processing
        # "--fast5_out", #!# if commented, does not save basecalled fast5 files
        "--progress_stats_frequency 1800" ## every 30 minutes
    shell:
        "guppy_basecaller --num_callers {threads} {params} "
        "--save_path {wildcards.tmp}01_processed_data/01_basecalling/{wildcards.run} "
        "--input_path {input} > {log} 2>&1; "
        "mkdir -p {wildcards.tmp}01_processed_data/01_basecalling/{wildcards.run}/guppy_basecaller_logs; "
        "mv {wildcards.tmp}01_processed_data/01_basecalling/{wildcards.run}/guppy_basecaller_log-* {wildcards.tmp}01_processed_data/01_basecalling/{wildcards.run}/guppy_basecaller_logs/"
