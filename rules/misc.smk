## file conversions

rule sam_bam:
    input:
        "01_processeddata/{run}/align/alignment.sam"
    output:
        "01_processeddata/{run}/align/alignment_sorted.bam"
    log:
        "01_processeddata/{run}/align/MeBaPiNa_sam_bam.log"
    benchmark:
        "01_processeddata/{run}/align/MeBaPiNa_sam_bam.benchmark.tsv"
    params:
        "-u" ## uncompressed bam to pipe
    conda:
        "../envs/samtools.yml"
    threads:
        config["machine"]["cpu"]
    shell:
        "samtools view --threads {threads} {params} {input} | samtools sort --threads {threads} -o {output} > {log} 2>&1"   
