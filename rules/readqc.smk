#############
## READ QC ##
#############

## TRIMM DEMULTIPLEX ##
#######################

rule trim_basecalled:
    ## input dependency is handled by the all_report rule
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
        "--epi2me", ## Use EPI2ME's demultiplexing algorithm (default: true)
        # "--guppy", ## Use Guppy's demultiplexing algorithm (default: false) #!# Multi threading only works with in guppy mode
        ("--kit RAB204" if BAC_KIT == "SQK-RAB204" else "--kit Auto") #!# otherwise there are no barcoding folders
    shell:
        "input_folder={wildcards.tmp}01_processed_data/01_basecalling/{wildcards.run}/pass/{wildcards.barc}; " ## input dependency is handled by the all rule
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
