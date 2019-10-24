    


rule minimap2:
    input:
        target="00_rawdata/reference_sequences/lambda/NC_001416.1_13-AUG-2018.mmi", # can be either genome index or genome fasta
        query="01_processeddata/{run}/pass"
    output:
        "01_processeddata/{run}/alignment/alignment.sam"
    log:
        "01_processeddata/{run}/alignment/MeBaPiNa_minimap2.log"
    benchmark:
        "01_processeddata/{run}/alignment/MeBaPiNa_minimap2.benchmark.tsv"
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
        "00_rawdata/reference_sequences/lambda/NC_001416.1_13-AUG-2018.fa"
    output:
        "00_rawdata/reference_sequences/lambda/NC_001416.1_13-AUG-2018.mmi"
    log:
        "00_rawdata/reference_sequences/lambda/MeBaPiNa_minimap2_index.log"
    benchmark:
        "00_rawdata/reference_sequences/lambda/MeBaPiNa_minimap2_index.benchmark.tsv"
    params:
        "-x map-ont" ## naopore specific
    conda:
        "../envs/minimap2.yml"
    threads:
        config["machine"]["cpu"]
    shell:
        "minimap2 -t {threads} {params} -d {output} {input} > {log} 2>&1"
