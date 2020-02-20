###############
## ALIGNMENT ##
###############

## ALIGNMENT ##
###############

rule splitting_filtered: #!# only required because of low memory avaialility
    input:
        "{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/filtered.fastq"
    output:
        temp(directory("{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/split"))
    log:
        "{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/MeBaPiNa_filtering.log"
    benchmark:
        "{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/MeBaPiNa_filtering.benchmark.tsv"
    shell:
        "mkdir {output}; "
        "awk 'BEGIN{{ file_num=0 }}; " ## split output into files with <=4000 sequences each
        "NR%16000==1{{ file_num++ }}; "
        "{{ print $0 > \"{output}/\" file_num \".fastq\" }}' {input} "
        ">> {log} 2>&1"

rule aligning_filtered:
    input:
        barc_dir="{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/split", 
        target=expand("{tmp}METADATA/Reference_Sequences/{reference}/reference.mmi", tmp = config["experiments"]["tmp"], reference = config["reference"]["source"])
    output:
        temp("{tmp}01_processed_data/03_alignment/{run}/{barc}/filtered.sam")
    log:
        "{tmp}01_processed_data/03_alignment/{run}/{barc}/MeBaPiNa_alignment.log"
    benchmark:
        "{tmp}01_processed_data/03_alignment/{run}/{barc}/MeBaPiNa_alignment.benchmark.tsv"
    params:
        "-x map-ont", ## naopore specific
        "-a" ## possition accurate CIGAR alignment in SAM output; much slower <- maybe skip?
        # -p FLOAT     min secondary-to-primary score ratio [0.8]
        # -N INT       retain at most INT secondary alignments [5]
    conda:
        "../envs/minimap2.yml"
    threads:
        8
    shell:
        "minimap2 -t {threads} {params} -o {output} "
        "{input.target} "
        "$(find {input.barc_dir} -type f -name \"*.fastq\") " #!# because of split input
        "> {log} 2>&1"

## FILE CONVERSION ##
#####################

rule converting_sam2bam:
    input:
        "{tmp}01_processed_data/03_alignment/{run}/{barc}/{altype}.sam"
    output:
        bam="{tmp}01_processed_data/03_alignment/{run}/{barc}/{altype}_sorted.bam",
        bai="{tmp}01_processed_data/03_alignment/{run}/{barc}/{altype}_sorted.bam.bai"
    log:
        "{tmp}01_processed_data/03_alignment/{run}/{barc}/MeBaPiNa_converting_{altype}.log"
    benchmark:
        "{tmp}01_processed_data/03_alignment/{run}/{barc}/MeBaPiNa_converting_{altype}.benchmark.tsv"
    conda:
        "../envs/samtools.yml"
    threads:
        2
    shell:
        "samtools sort --threads {threads} -o {output.bam} {input} > {log} 2>&1; "
        "samtools index -@ {threads} {output.bam} >> {log} 2>&1"

## CALIBRATION STRAND ##
########################

rule aligning_calibration_strands:
    input:
        calib_dir="{tmp}01_processed_data/01_basecalling/{run}/calibration_strands", 
        target="{tmp}METADATA/Reference_Sequences/lambda_3.6kb/reference.mmi"
    output:
        temp("{tmp}01_processed_data/03_alignment/{run}/lambda/calibration.sam")
    log:
        "{tmp}01_processed_data/03_alignment/{run}/lambda/MeBaPiNa_alignment.log"
    benchmark:
        "{tmp}01_processed_data/03_alignment/{run}/lambda/MeBaPiNa_alignment.benchmark.tsv"
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
