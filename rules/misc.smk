
rule sam_to_bam:
    input:
        "01_processeddata/{run}/{align}/{barc}_alignment.sam"
    output:
        bam="01_processeddata/{run}/{align}/{barc}_alignment_sorted.bam",
        bai="01_processeddata/{run}/{align}/{barc}_alignment_sorted.bam.bai"
    log:
        "01_processeddata/{run}/{align}/{barc}_MeBaPiNa_sam_bam.log"
    benchmark:
        "01_processeddata/{run}/{align}/{barc}_MeBaPiNa_sam_bam.benchmark.tsv"
    version:
        subprocess.check_output("samtools --version | awk 'NR==1{print $NF}'", shell=True)
    params:
        "-u" ## uncompressed bam to pipe
    conda:
        "../envs/samtools.yml"
    threads:
        config["machine"]["cpu"]
    shell:
        "samtools view --threads {threads} {params} {input} | "
        "samtools sort --threads {threads} -o {output.bam} > {log} 2>&1; "
        "samtools index -@ {threads} {output.bam} > {log} 2>&1"

rule fasq_pipe:
    input:
        "01_processeddata/{run}/basecall/pass"
    output:
        temp("02_analysis/{run}/basecall/nanoqc/pipe.fastq.gz") ## pipe didn't work
    shell:
        "cat $(find {input} -type f -name \"*.fastq.gz\") "
        ">> {output}"

rule find_reads_in_fastq:
    input:
        "01_processeddata/{run}/basecall/pass"
    output:
        "01_processeddata/{run}/basecall/find_reads_in_fastq.txt"
    log:
        "01_processeddata/{run}/basecall/MeBaPiNa_find_fastq_in_fast5.log"
    benchmark:
        "01_processeddata/{run}/basecall/MeBaPiNa_find_fastq_in_fast5.benchmark.tsv"
    conda:
        "../envs/find_fastq_in_fast5.yml"
    script:
        "../scripts/find_fastq_in_fast5.py"

rule split_seqsum_barc:
    input:
        "01_processeddata/{run}/basecall/sequencing_summary.txt"
    output:
        directory("01_processeddata/{run}/basecall/sequencing_summary_barc")
    log:
        "01_processeddata/{run}/basecall/sequencing_summary_barc/MeBaPiNa_split_seqsum_barc.log"
    benchmark:
        "01_processeddata/{run}/basecall/sequencing_summary_barc/MeBaPiNa_split_seqsum_barc.benchmark.tsv"
    conda:
        "../envs/pycoqc.yml"
    params:
        "--output_unclassified", ## keep unclassified in sparate file
        "--min_barcode_percent 0.001", ## barcodes below 0.001% of reads are removed (1 in 100'000)
        "--verbose"
    shell:
        "Barcode_split {params} "
        "--summary_file {input} "
        "--output_dir {output} > {log} 2>&1"

def input_aggregate(wildcards):
    from os import listdir
    basecall_dir = checkpoints.guppy.get(run=wildcards.run).output[0]
    barcs = listdir(basecall_dir)
    barc_dirs = expand("02_analysis/{run}/align/{barc}_pycoqc/pycoQC_report.json",
        run=wildcards.run,
        barc=barcs)
    return barc_dirs

rule aggr_align_barc:
    input:
        input_aggregate
    output:
        "02_analysis/{run}/align/MeBaPiNa_barcode_aggregation.txt"
    shell:
        "echo {input} > {output}"
