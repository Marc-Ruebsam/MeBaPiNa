    


rule minimap2:
    input:
        target="01_processeddata/reference_sequences/lambda/some.mmi", # can be either genome index or genome fasta
        query="01_processeddata/{run}/pass/{reads}.fastq.gz"
    output:
        "01_processeddata/{run}/alignment/lambda_reference.sam"
    log:
        "01_processeddata/{run}/alignment/MeBaPiNa_minimap2.log"
    benchmark:
        "01_processeddata/{run}/alignment/MeBaPiNa_minimap2.benchmark.tsv"
    params:
        "-x map-ont" ## naopore specific
        "-a" ## possition accurate CIGAR alignment in SAM output; much slower <- maybe skip?
    threads:
        config["machine"]["cpu"]
    shell:
        "minimap2 -t {threads} {params} -o {output} {input.target} {input.query} > {log} 2>&1"

rule minimap2_index:
    input:
        "01_processeddata/reference_sequences/lambda/some.fa"
    output:
        "01_processeddata/reference_sequences/lambda/some.mmi"
    log:
        "01_processeddata/reference_sequences/lambda/MeBaPiNa_minimap2_index.log"
    benchmark:
        "01_processeddata/reference_sequences/lambda/MeBaPiNa_minimap2_index.benchmark.tsv"
    params:
        "-x map-ont" ## naopore specific
    threads:
        config["machine"]["cpu"]
    shell:
        "minimap2 -t {threads} {params} -d {output} {input} > {log} 2>&1"
