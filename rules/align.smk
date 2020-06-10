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
        "mkdir {output[0]}; "
        "touch {log}"

rule aligning_filtered:
    input:
        barc_dir="{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/split",
        target="{tmp}METADATA/Reference_Sequences/{reference}/reference.mmi"
    output:
        temp("{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/aligned.sam")
    log:
        "{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/MeBaPiNa_aligning_filtered.log"
    benchmark:
        "{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/MeBaPiNa_aligning_filtered.benchmark.tsv"
    shell:
        "touch {output[0]}; "
        "lg={log}; "
        "cat ${{lg/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} > {log}"

## FILE CONVERSION ##
#####################

rule filter_aligned:
    input:
        "{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/aligned.sam"
    output:
        sam=temp("{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/filtered.sam"),
        bam="{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/filteredsorted.bam",
        bai="{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/filteredsorted.bam.bai"
    log:
        "{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/MeBaPiNa_filter_aligned.log"
    benchmark:
        "{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/MeBaPiNa_filter_aligned.benchmark.tsv"
    shell:
        "touch {output[0]}; "
        "out2={output[1]}; "
        "cp ${{out2/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} ${{out2}}; "
        "out3={output[2]}; "
        "cp ${{out3/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} ${{out3}}; "
        "lg={log}; "
        "cat ${{lg/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} > {log}"

rule convert_bambacktosam: #!# used in case of rerunning
    input:
        bam="{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/filteredsorted.bam"
    output:
        sam=temp("{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/filtered.sam")
    log:
        "{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/MeBaPiNa_bambacktosam.log"
    benchmark:
        "{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/MeBaPiNa_bambacktosam.benchmark.tsv"
    shell:
        "touch {output[0]}; "
        "lg={log}; "
        "cat ${{lg/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} > {log}"
ruleorder: convert_bambacktosam > filter_aligned

rule counttax_aligned:
    input:
        sam="{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/filtered.sam",
        kronataxlist="{tmp}METADATA/Reference_Sequences/{reference}/krona/{reftype}/taxlist.txt",
        kronaseq2tax="{tmp}METADATA/Reference_Sequences/{reference}/krona/{reftype}/seqid2taxid.map"
    output:
        counttaxlist="{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}_{reftype}/aligned.counttaxlist"
    log:
        "{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_counttax_aligned.log"
    benchmark:
        "{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_counttax_aligned.benchmark.tsv"
    shell:
        "out1={output[0]}; "
        "cp ${{out1/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} ${{out1}}; "
        "lg={log}; "
        "cat ${{lg/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} > {log}"
