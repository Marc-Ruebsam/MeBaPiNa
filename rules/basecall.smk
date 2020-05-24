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
        list(filter(None,[
        "{tmp}01_processed_data/01_basecalling/{run}/sequencing_summary.txt",
        directory("{tmp}01_processed_data/01_basecalling/{run}/pass"),
        (directory("{tmp}01_processed_data/01_basecalling/{run}/calibration_strands") if LAM_DCS else "") ## evaluation of Lambda calibration strands only when specified
        ]))
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
        # "--compress_fastq",
        "--fast5_out", #!#
        "--progress_stats_frequency 1800" ## every 30 minutes
    shell:
        "guppy_basecaller --num_callers {threads} {params} "
        "--save_path {wildcards.tmp}01_processed_data/01_basecalling/{wildcards.run} "
        "--input_path {input} > {log} 2>&1; "
        "mkdir -p {wildcards.tmp}01_processed_data/01_basecalling/{wildcards.run}/guppy_basecaller_logs; "
        "mv {wildcards.tmp}01_processed_data/01_basecalling/{wildcards.run}/guppy_basecaller_log-* {wildcards.tmp}01_processed_data/01_basecalling/{wildcards.run}/guppy_basecaller_logs/"

## TRIMM DEMULTIPLEX ##
#######################

rule trim_basecalled:
    ## input dependency is handled by the all rule
    output:
        "{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/trimmed.fastq"
    log:
        "{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/MeBaPiNa_trim_basecalled.log"
    benchmark:
        "{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/MeBaPiNa_trim_basecalled.benchmark.tsv"
    conda:
        "../envs/qcat.yml"
    threads:
        1
    params:
        "--min-score 70", ## Minimum barcode score. Barcode calls with a lower score will be discarded. Must be between 0 and 100. (default: 60)
        "--detect-middle", ## Search for adapters in the whole read #!# slower
        # "--min-read-length " + config["filtering"]["len_min"], ## Reads short than <min-read-length> after trimming will be discarded. #!# Should be done in next step
        "--trim", ## Remove adapter and barcode sequences from reads.
        ("--kit RAB204" if BAC_KIT == "SQK-RAB204" else "--kit Auto") #!# otherwise there are no barcoding folders
    shell:
        "input_folder={wildcards.tmp}01_processed_data/01_basecalling/{wildcards.run}/pass/{barc}; " ## input dependency is handled by the all rule
        "output_folder={wildcards.tmp}01_processed_data/02_trimming_filtering/{wildcards.run}/{wildcards.barc}; " ## directory name of barcode
        "find ${{input_folder}} -type f -name \"*.fastq\" -exec cat {{}} \\; | " ## paste whole read content of all fast files to std-in of qcat
        "qcat --threads {threads} {params} "
        "--barcode_dir ${{output_folder}} "
        "> {log} 2>&1; "
        "mv ${{output_folder}}/{wildcards.barc}.fastq {output} >> {log} 2>&1; " ## rename barcode fastq to "trimmed.fastq"
        "mkdir -p ${{output_folder}}/others >> {log} 2>&1; " ## create folder for all other barcode files sorted out during demultiplexing
        "find ${{output_folder}} -type f \( -name \"*barcode*\" -o -name \"*none*\" \) " ## move other barcodes to "others" directory
        "-exec mv {{}} ${{output_folder}}/others \; >> {log} 2>&1"

## FILTER ##
############

rule filter_trimmed:
    input:
        fastq="{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/trimmed.fastq",
        seqsum="{tmp}01_processed_data/01_basecalling/{run}/sequencing_summary.txt"
    output:
        "{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/filtered.fastq"
    log:
        "{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/MeBaPiNa_filter_trimmed.log"
    benchmark:
        "{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/MeBaPiNa_filter_trimmed.benchmark.tsv"
    conda:
        "../envs/nanopack.yml"
    params:
        ("--length " + config["filtering"]["len_min"]),
        ("--maxlength " + config["filtering"]["len_max"]),
        ("--quality " + config["filtering"]["q_min"])
    shell:
        "NanoFilt {params} --summary {input.seqsum} "
        "--logfile {log} {input.fastq}"
        "> {output} 2> {log}"
