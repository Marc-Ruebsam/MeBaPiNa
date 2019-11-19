
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
        "01_processeddata/{run}/basecall/sequencing_summary.txt" ## used as dummy for the other folders
    output:
        temp("02_analysis/{run}/basecall/nanoqc/pipe.fastq.gz") ## pipe didn't work
    shell:
        "cat $(find 01_processeddata/{wildcards.run}/basecall/pass -type f -name \"*.fastq.gz\") "
        ">> {output}"
