#################
## BASECALLING ##
#################

## HOUSE KEEPING ##
###################

rule moving_raw:
    output:
        directory("{tmp}00_raw_data/{run}/fast5")
    shell:
        "long_path=$(find {wildcards.tmp}00_raw_data/ -name \"{wildcards.run}\"); "
        "mv ${{long_path}} \"{wildcards.tmp}00_raw_data/{wildcards.run}\"; "
        "find {wildcards.tmp}00_raw_data/ -depth -type d -empty -delete" #!# delete empty directories

rule report_moving_raw:
    input:
        "{tmp}00_raw_data/{run}/fast5"
    output:
        temp("{tmp}00_raw_data/{run}/MeBaPiNa_moving_raw.report")
    params:
        " ".join(SAMPLES.values())
    shell:
        "all_IDs=( $(echo \"{params}\") ); "
        "for id in ${{all_IDs[@]}}; do echo "
        ## "Sample name;File/directory;Completion date;Checksum;Performed by;Description"
        "\"${{id}};{input};$(stat -c %y {input});NA;MeBaPiNa;Moved raw fast5 files to path ready for analysis.\" "
        ">> {output}; done"

rule report_basecalling_raw:
    input:
        list(filter(None,[
        "{tmp}01_processed_data/01_basecalling/{run}/sequencing_summary.txt",
        directory("{tmp}01_processed_data/01_basecalling/{run}/pass"),
        (directory("{tmp}01_processed_data/01_basecalling/{run}/calibration_strands") if LAM_DCS else "") ## evaluation of Lambda calibration strands only when specified
        ]))
    output:
        temp("{tmp}00_raw_data/{run}/MeBaPiNa_basecalling_raw.report")
    params:
        " ".join(SAMPLES.values())
    shell:
        "all_IDs=( $(echo \"{params}\") ); "
        "for id in ${{all_IDs[@]}}; do echo "
        ## "Sample name;File/directory;Completion date;Checksum;Performed by;Description"
        "\"${{id}};{input};$(stat -c %y {input[0]});NA;MeBaPiNa;Basecalled and demultiplexed raw fast5 files into fastq.\" "
        ">> {output}; done"

rule report_trimming_basecalled:
    input:
        "{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/trimmed.fastq"
    output:
        temp("{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/MeBaPiNa_trimming_basecalled.report")
    shell:
        "echo "
        ## "Sample name;File/directory;Completion date;Checksum;Performed by;Description"
        "\"{wildcards.barc};{input};$(stat -c %y {input});NA;MeBaPiNa;Trimmed adapters and barcodes from reads and demultiplexed a second time.\" "
        ">> {output}"

rule report_filtering_trimmed:
    input:
        "{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/filtered.fastq"
    output:
        temp("{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/MeBaPiNa_filtering_trimmed.report")
    shell:
        "echo "
        ## "Sample name;File/directory;Completion date;Checksum;Performed by;Description"
        "\"{wildcards.barc};{input};$(stat -c %y {input});NA;MeBaPiNa;Length and quality filtered reads.\" "
        ">> {output}"


## BASECALL DEMULTIPLEX ##
##########################

checkpoint basecalling_raw:
    input:
        "{tmp}00_raw_data/{run}/fast5"
    output:
        list(filter(None,[
        "{tmp}01_processed_data/01_basecalling/{run}/sequencing_summary.txt",
        directory("{tmp}01_processed_data/01_basecalling/{run}/pass"),
        (directory("{tmp}01_processed_data/01_basecalling/{run}/calibration_strands") if LAM_DCS else "") ## evaluation of Lambda calibration strands only when specified
        ]))
    log:
        "{tmp}01_processed_data/01_basecalling/{run}/MeBaPiNa_basecalling_raw.log"
    benchmark:
        "{tmp}01_processed_data/01_basecalling/{run}/MeBaPiNa_basecalling_raw.benchmark.tsv"
    version:
        subprocess.check_output("guppy_basecaller --version | awk '{print $NF}'", shell=True)
    threads:
        22
    params:
        "--flowcell " + FLOWCELL,
        "--kit " + SEQ_KIT,
        "--barcode_kits " + BAC_KIT, ## always includes demultiplexing (all reads marked as unclassified if no barcodes were used)
        ("--calib_detect" if LAM_DCS else ""), ## includes detection of lambda clibration strands
        # "--min_score 75", ## minimum score for barcode match (default 60)
        # # "--require_barcodes_both_ends", ## very stringent
        # # "--detect_mid_strand_barcodes", ## very slow
        # # "--min_score_mid_barcodes 70"
        # # "--min_score_rear_override 70"
        # # "--trim_barcodes",
        # # " --num_extra_bases_trim 10",
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

rule trimming_basecalled:
    ## input dependency is handled by the all rule
    output:
        "{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/trimmed.fastq"
    log:
        "{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/MeBaPiNa_trimming_basecalled.log"
    benchmark:
        "{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/MeBaPiNa_trimming_basecalled.benchmark.tsv"
    conda:
        "../envs/qcat.yml"
    threads:
        1
    params:
        "--min-score 70", ## Minimum barcode score. Barcode calls with a lower score will be discarded. Must be between 0 and 100. (default: 60)
        "--detect-middle", #!# Search for adapters in the whole read
        # "--min-read-length " + config["filtering"]["len_min"], ## Reads short than <min-read-length> after trimming will be discarded.
        "--trim", ## Remove adapter and barcode sequences from reads.
        ("--kit RAB204" if BAC_KIT == "SQK-RAB204" else "--kit Auto") #!#
    shell:
        "input_folder={wildcards.tmp}01_processed_data/01_basecalling/{wildcards.run}/pass; " ## input dependency is handled by the all rule
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

rule filtering_trimmed:
    input:
        fastq="{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/trimmed.fastq",
        seqsum="{tmp}01_processed_data/01_basecalling/{run}/sequencing_summary.txt"
    output:
        "{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/filtered.fastq"
    log:
        "{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/MeBaPiNa_filtering_trimmed.log"
    benchmark:
        "{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/MeBaPiNa_filtering_trimmed.benchmark.tsv"
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
