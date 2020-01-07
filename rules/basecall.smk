##########################
## BASECALL DEMULTIPLEX ##
##########################

checkpoint guppy:
    input:
        "/mnt/NRD/{run}/fast5"
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
        22
    params:
        ("--flowcell " + config["guppy"]["flowcell"]),
        ("--kit " + config["guppy"]["seq_kit"]),
        ("--calib_detect" if config["guppy"]["lam_DCS"] else ""), ## includes detection of lambda clibration strands
        "--barcode_kits " + BAC_KIT, ## always includes demultiplexing (all reads marked as unclassified if no barcodes were used)
        # # "--require_barcodes_both_ends",
        # "--detect_mid_strand_barcodes",
        # "--min_score 75",
        # "--min_score_rear_override 70"
        # "--min_score_mid_barcodes 70"
        # "--trim_barcodes",
        # " --num_extra_bases_trim 10",
        "--qscore_filtering",
        ("--min_qscore " + config["guppy"]["q_min"]),
        ("--device cuda:all:100% "
        + "--gpu_runners_per_device 6 " 
        + "--chunks_per_runner 1536 " 
        + "--chunk_size 1000 " 
        + "--chunks_per_caller 10000 " 
        + "--num_barcode_threads 8" 
        if config["machine"]["gpu"] else ""),
        # "--compress_fastq",
        "--fast5_out",
        "--progress_stats_frequency 1800"
    shell:
        "guppy_basecaller --num_callers {threads} {params} "
        "--save_path 01_processeddata/{wildcards.run}/basecall "
        "--input_path {input} > {log} 2>&1; "
        "mkdir -p 01_processeddata/{wildcards.run}/basecall/guppy_basecaller_logs; "
        "mv 01_processeddata/{wildcards.run}/basecall/guppy_basecaller_log-* 01_processeddata/{wildcards.run}/basecall/guppy_basecaller_logs/"

#######################
## TRIMM DEMULTIPLEX ##
#######################

rule qcat:
    input:
        "01_processeddata/{run}/basecall/pass/{barc}"
    output:
        directory("01_processeddata/{run}/trim/{barc}")
    log:
        "01_processeddata/{run}/trim/{barc}_MeBaPiNa_qcat.log"
    benchmark:
        "01_processeddata/{run}/trim/{barc}_MeBaPiNa_qcat.benchmark.tsv"
    conda:
        "../envs/qcat.yml"
    threads:
        1
    params:
        "--min-score 70", ## Minimum barcode score. Barcode calls with a lower score will be discarded. Must be between 0 and 100. (default: 60)
        "--detect-middle", ## Search for adapters in the whole read
        ("--min-read-length " + config["guppy"]["len_min"]), ## Reads short than <min-read-length> after trimming will be discarded.
        "--trim", ## Remove adapter and barcode sequences from reads.
        ("--kit RAB204" if BAC_KIT == "SQK-RAB204" else "--kit Auto")
    shell:
        "find {input} -type f -name \"*.fastq\" -exec cat {{}} \\; | "
        "qcat --threads {threads} {params} "
        "--barcode_dir {output} "
        "> {log} 2>&1"

############
## FILTER ##
############

rule nanofilt:
    input:
        fastq="01_processeddata/{run}/trim/{barc}",
        seqsum="01_processeddata/{run}/basecall/sequencing_summary.txt"
    output:
        "01_processeddata/{run}/filter/{barc}.fastq"
    log:
        "01_processeddata/{run}/filter/{barc}_MeBaPiNa_nanofilt.log"
    benchmark:
        "01_processeddata/{run}/filter/{barc}_MeBaPiNa_nanofilt.benchmark.tsv"
    conda:
        "../envs/nanopack.yml"
    params:
        ("--length " + config["guppy"]["len_min"]),
        ("--maxlength " + config["guppy"]["len_max"]),
        ("--quality " + config["guppy"]["q_min"])
    shell:
        "NanoFilt {params} --summary {input.seqsum} "
        "--logfile {log} {input.fastq}/{wildcards.barc}.fastq"
        "> {output} 2> {log}"
