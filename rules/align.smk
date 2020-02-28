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
        temp("{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/align.sam")
    log:
        "{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/MeBaPiNa_alignment.log"
    benchmark:
        "{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/MeBaPiNa_alignment.benchmark.tsv"
    conda:
        "../envs/minimap2.yml"
    params:
        "-x map-ont", ## naopore specific
        "-a", ## possition accurate CIGAR alignment in SAM output; much slower <- maybe use -c for PAF instead?
        "-p 1", ## only retain multi mappings with same highest score
        "-N 1" ## one secondary alignment is enough to identify multimapping reads
    threads:
        34
    shell:
        "minimap2 -t {threads} {params} -o {output} "
        "{input.target} "
        "$(find {input.barc_dir} -type f -name \"*.fastq\") " #!# because of split input
        "> {log} 2>&1"

## FILE CONVERSION ##
#####################

rule filter_aligned:
    input:
        "{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/{altype}.sam"
    output:
        sam="{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/{altype}_filtered.sam",
        bam="{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/{altype}_filteredsorted.bam",
        bai="{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/{altype}_filteredsorted.bam.bai"
    log:
        "{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/MeBaPiNa_converting_{altype}.log"
    benchmark:
        "{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/MeBaPiNa_converting_{altype}.benchmark.tsv"
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
        temp("{tmp}01_processed_data/03_alignment/{run}/{reference}/lambda/calibration.sam")
    log:
        "{tmp}01_processed_data/03_alignment/{run}/{reference}/lambda/MeBaPiNa_alignment.log"
    benchmark:
        "{tmp}01_processed_data/03_alignment/{run}/{reference}/lambda/MeBaPiNa_alignment.benchmark.tsv"
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
