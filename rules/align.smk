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
        "{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/MeBaPiNa_splitting_filtered.log"
    benchmark:
        "{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/MeBaPiNa_splitting_filtered.benchmark.tsv"
    shell:
        "mkdir -p {output}; "
        "awk 'BEGIN{{ file_num=0 }}; " ## split output into files with <=4000 sequences each
        "NR%16000==1{{ file_num++ }}; "
        "{{ print $0 > \"{output}/\" file_num \".fastq\" }}' {input} "
        ">> {log} 2>&1"

rule aligning_filtered:
    input:
        barc_dir="{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/split",
        target=expand("{tmp}METADATA/Reference_Sequences/{reference}/reference.mmi", tmp = config["experiments"]["tmp"], reference = config["reference"]["source"])
    output:
        temp("{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/aligned.sam")
    log:
        "{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/MeBaPiNa_aligning_filtered.log"
    benchmark:
        "{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/MeBaPiNa_aligning_filtered.benchmark.tsv"
    conda:
        "../envs/minimap2.yml"
    params:
        "-x map-ont", ## naopore specific
        "-a", ## possition accurate CIGAR alignment in SAM output; much slower <- maybe use -c for PAF instead?
        "-p 1", ## only retain multi mappings with same highest score
        "-N 1" ## one secondary alignment is enough to identify multimapping reads
    threads:
        8
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
        sam=temp("{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/{altype}_filtered.sam"),
        bam="{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/{altype}_filteredsorted.bam",
        bai="{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/{altype}_filteredsorted.bam.bai"
    log:
        "{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/MeBaPiNa_filter_{altype}.log"
    benchmark:
        "{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/MeBaPiNa_filter_{altype}.benchmark.tsv"
    conda:
        "../envs/samtools.yml"
    params:
        "-F 2048", ## exclude chimeric/supplementary alignments
        "-q 1" ## MAPQ >= 1 filters out all multimapping reads as they always have MAPQ=0
    threads:
        4
    shell:
        "samtools view -h -@ $(({threads} / 2)) {params} 2> {log} {input} | "
        "tee {output.sam} | "
        "samtools sort --threads $(({threads} / 2)) -o {output.bam} >> {log} 2>&1; "
        "samtools index -@ {threads} {output.bam} >> {log} 2>&1"

rule convert_bambacktosam:
    input:
        bam="{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/{altype}_filteredsorted.bam"
    output:
        sam="{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/{altype}_filtered.sam"
    log:
        "{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/MeBaPiNa_filter_{altype}.log"
    benchmark:
        "{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/MeBaPiNa_filter_{altype}.benchmark.tsv"
    conda:
        "../envs/samtools.yml"
    threads:
        1
    shell:
        "samtools view -h -@ $(({threads} / 2)) {params} -o {output} {input} > {log} 2>&1"
ruleorder: convert_bambacktosam > filter_aligned

rule counttax_aligned:
    input:
        sam="{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/aligned_filtered.sam",
        kronataxlist="{tmp}METADATA/Reference_Sequences/{reference}/krona/{reftype}/taxlist.txt",
        kronaseq2tax="{tmp}METADATA/Reference_Sequences/{reference}/krona/{reftype}/seqid2taxid.map"
    output:
        counttaxlist="{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}_{reftype}/aligned.counttaxlist"
    log:
        "{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_counttax_aligned.log"
    benchmark:
        "{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_counttax_aligned.benchmark.tsv"
    conda:
        "../envs/python.yml"
    script:
        "../scripts/convert_sam.py" ## filters by read length and average base diversion

## CALIBRATION STRAND ##
########################

rule aligning_calibstr:
    input:
        calib_dir="{tmp}01_processed_data/01_basecalling/{run}/calibration_strands",
        target="{tmp}METADATA/Reference_Sequences/lambda_3.6kb/reference.mmi"
    output:
        temp("{tmp}01_processed_data/03_alignment/{run}/{reference}/lambda/calibration.sam")
    log:
        "{tmp}01_processed_data/03_alignment/{run}/{reference}/lambda/MeBaPiNa_aligning_calibstr.log"
    benchmark:
        "{tmp}01_processed_data/03_alignment/{run}/{reference}/lambda/MeBaPiNa_aligning_calibstr.benchmark.tsv"
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
