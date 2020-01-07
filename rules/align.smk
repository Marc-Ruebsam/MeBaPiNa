rule minimap2_index:
    input:
        "00_rawdata/reference_sequences/{reference}.fasta"
    output:
        "00_rawdata/reference_sequences/{reference}.mmi"
    log:
        "00_rawdata/reference_sequences/{reference}_MeBaPiNa_minimap2_index.log"
    benchmark:
        "00_rawdata/reference_sequences/{reference}_MeBaPiNa_minimap2_index.benchmark.tsv"
    params:
        "-x map-ont" ## naopore specific
    conda:
        "../envs/minimap2.yml"
    threads:
        2
    shell:
        "minimap2 -t {threads} {params} -d {output} {input} > {log} 2>&1"

rule fastq_split_minimap2:
    input:
        "01_processeddata/{run}/filter/{barc}.fastq"
    output:
        temp(directory("01_processeddata/{run}/filter/{barc}_split"))
    log:
        "01_processeddata/{run}/filter/{barc}_MeBaPiNa_split.log"
    benchmark:
        "01_processeddata/{run}/filter/{barc}_MeBaPiNa_split.tsv"
    shell:
        "mkdir {output}; "
        "awk 'BEGIN{{ file_num=0 }}; " ## split output into files with 4000 sequences each
        "NR%16000==1{{ file_num++ }}; "
        "{{ print $0 > \"{output}/{wildcards.barc}_\" file_num \".fastq\" }}' {input} "
        "> {log} 2>&1"

rule minimap2:
    input:
        barc_dir="01_processeddata/{run}/filter/{barc}_split", 
        target=expand("00_rawdata/reference_sequences/{reference}.mmi", reference=config["align"]["reference"])
    output:
        temp("01_processeddata/{run}/align/{barc}_alignment.sam")
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
        8
    shell:
        "minimap2 -t {threads} {params} -o {output} "
        "{input.target} "
        "$(find {input.barc_dir} -type f -name \"*.fastq\") "
        "> {log} 2>&1"

rule minimap2_calibration_strands:
    input:
        calib_dir="01_processeddata/{run}/basecall/calibration_strands", 
        target="00_rawdata/reference_sequences/lambda/lambda_3.6kb.mmi"
    output:
        "01_processeddata/{run}/align_calibration_strands/lambda_alignment.sam"
    log:
        "01_processeddata/{run}/align_calibration_strands/lambda_MeBaPiNa_minimap2.log"
    benchmark:
        "01_processeddata/{run}/align_calibration_strands/lambda_MeBaPiNa_minimap2.benchmark.tsv"
    params:
        "-x map-ont", ## naopore specific
        "-a" ## possition accurate CIGAR alignment in SAM output; much slower <- maybe skip?
    conda:
        "../envs/minimap2.yml"
    threads:
        8
    shell:
        "minimap2 -t {threads} {params} -o {output} "
        "{input.target} "
        "$(find {input.calib_dir} -type f -name \"*.fastq\") "
        "> {log} 2>&1"
