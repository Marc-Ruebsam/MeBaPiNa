rule minimap2:
    input:
        target="00_rawdata/reference_sequences/lambda/{reference}.fasta",
        query="01_processeddata/{run}/basecall/pass"
    output:
        "01_processeddata/{run}/align/alignment_{reference}.sam"
    log:
        "01_processeddata/{run}/align/alignment_{reference}_MeBaPiNa_minimap2.log"
    benchmark:
        "01_processeddata/{run}/align/alignment_{reference}_MeBaPiNa_minimap2.benchmark.tsv"
    params:
        "-x map-ont", ## naopore specific
        "-a" ## possition accurate CIGAR alignment in SAM output; much slower <- maybe skip?
    conda:
        "../envs/minimap2.yml"
    threads:
        config["machine"]["cpu"]
    shell:
        "minimap2 -t {threads} {params} -o {output} {input.target} {input.query}/* > {log} 2>&1"

rule minimap2_index:
    input:
        "00_rawdata/reference_sequences/lambda/{reference}.fasta"
    output:
        "00_rawdata/reference_sequences/lambda/{reference}.mmi"
    log:
        "00_rawdata/reference_sequences/lambda/{reference}_MeBaPiNa_minimap2_index.log"
    benchmark:
        "00_rawdata/reference_sequences/lambda/{reference}_MeBaPiNa_minimap2_index.benchmark.tsv"
    params:
        "-x map-ont" ## naopore specific
    conda:
        "../envs/minimap2.yml"
    threads:
        config["machine"]["cpu"]
    shell:
        "minimap2 -t {threads} {params} -d {output} {input} > {log} 2>&1"

rule minimap2_lam:
    input:
        target="00_rawdata/reference_sequences/lambda/lambda_3.6kb.fasta",
        query="01_processeddata/{run}/calibration_strands"
    output:
        "01_processeddata/{run}/align_lam/alignment.sam"
    log:
        "01_processeddata/{run}/align_lam/MeBaPiNa_minimap2.log"
    benchmark:
        "01_processeddata/{run}/align_lam/MeBaPiNa_minimap2.benchmark.tsv"
    params:
        "-x map-ont", ## naopore specific
        "-a" ## possition accurate CIGAR alignment in SAM output; much slower <- maybe skip?
    conda:
        "../envs/minimap2.yml"
    threads:
        config["machine"]["cpu"]
    shell:
        "minimap2 -t {threads} {params} -o {output} {input.target} {input.query}/* > {log} 2>&1"
