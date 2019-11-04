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

rule minimap2:
    input:
        barc_dir="01_processeddata/{run}/basecall/pass/{barc}", 
        target=expand("00_rawdata/reference_sequences/lambda/{reference}.mmi", reference=config["align"]["reference"])
    output:
        "01_processeddata/{run}/align/{barc}_alignment.sam"
    log:
        "01_processeddata/{run}/align/{barc}_MeBaPiNa_minimap2.log"
    benchmark:
        "01_processeddata/{run}/align/{barc}_MeBaPiNa_minimap2.benchmark.tsv"
    params:
        "-x map-ont", ## naopore specific
        "-a" ## possition accurate CIGAR alignment in SAM output; much slower <- maybe skip?
    conda:
        "../envs/minimap2.yml"
    threads:
        config["machine"]["cpu"]
    shell:
        "minimap2 -t {threads} {params} -o {output} "
        "{input.target} "
        "$(find {input.barc_dir} -type f -name \"*.fastq.gz\") "
        "> {log} 2>&1"

rule minimap2_calibration_strands:
    input:
        calib_dir="01_processeddata/{run}/basecall/calibration_strands", 
        target="00_rawdata/reference_sequences/lambda/lambda_3.6kb.mmi"
    output:
        "01_processeddata/{run}/align_calibration_strands/alignment.sam"
    log:
        "01_processeddata/{run}/align_calibration_strands/MeBaPiNa_minimap2.log"
    benchmark:
        "01_processeddata/{run}/align_calibration_strands/MeBaPiNa_minimap2.benchmark.tsv"
    params:
        "-x map-ont", ## naopore specific
        "-a" ## possition accurate CIGAR alignment in SAM output; much slower <- maybe skip?
    conda:
        "../envs/minimap2.yml"
    threads:
        config["machine"]["cpu"]
    shell:
        "minimap2 -t {threads} {params} -o {output} "
        "{input.target} "
        "$(find {input.calib_dir} -type f -name \"*.fastq.gz\") "
        "> {log} 2>&1"
