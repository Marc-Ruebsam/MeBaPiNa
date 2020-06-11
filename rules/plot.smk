##############
## PLOTTING ##
##############

## BASECALLING ##
#################

## NANOPLOT ##

rule plot_nanoplot_seqsum_basecall:
    input:
        "{tmp}01_processed_data/01_basecalling/{run}/sequencing_summary.txt"
    output:
        "{tmp}02_analysis_results/01_basecalling/{run}/nanoplot/NanoStats.txt",
        "{tmp}02_analysis_results/01_basecalling/{run}/nanoplot/NanoPlot-report.html"
    log:
        "{tmp}02_analysis_results/01_basecalling/{run}/nanoplot/MeBaPiNa_nanoplot_seqsum_basecall.log"
    benchmark:
        "{tmp}02_analysis_results/01_basecalling/{run}/nanoplot/MeBaPiNa_nanoplot_seqsum_basecall.benchmark.tsv"
    shell:
        "out1={output[0]}; "
        "cp ${{out1/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} ${{out1}}; "
        "out2={output[1]}; "
        "cp ${{out2/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} ${{out2}}; "
        "lg={log}; "
        "cat ${{lg/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} > {log}"

## PYCOQC ##

rule plot_pycoqc_seqsum_basecall:
    input:
        "{tmp}01_processed_data/01_basecalling/{run}/sequencing_summary/sequencing_summary_downsampled.txt" #!# "{tmp}01_processed_data/01_basecalling/{run}/sequencing_summary/sequencing_summary_sorted.txt"
    output:
        html="{tmp}02_analysis_results/01_basecalling/{run}/pycoqc/pycoQC_report.html",
        json="{tmp}02_analysis_results/01_basecalling/{run}/pycoqc/pycoQC_report.json"
    log:
        "{tmp}02_analysis_results/01_basecalling/{run}/pycoqc/MeBaPiNa_pycoqc_seqsum_basecall.log"
    benchmark:
        "{tmp}02_analysis_results/01_basecalling/{run}/pycoqc/MeBaPiNa_pycoqc_seqsum_basecall.benchmark.tsv"
    shell:
        "out1={output[0]}; "
        "cp ${{out1/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} ${{out1}}; "
        "out2={output[1]}; "
        "cp ${{out2/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} ${{out2}}; "
        "lg={log}; "
        "cat ${{lg/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} > {log}"

## NANOCOMP ##

rule plot_nanocomp_seqsum_basecall:
    input:
        "{tmp}01_processed_data/01_basecalling/{run}/sequencing_summary/sequencing_summary_sorted.txt"
    output:
        "{tmp}02_analysis_results/01_basecalling/{run}/nanocomp/NanoStats.txt",
        "{tmp}02_analysis_results/01_basecalling/{run}/nanocomp/NanoComp-report.html"
    log:
        "{tmp}02_analysis_results/01_basecalling/{run}/nanocomp/MeBaPiNa_nanocomp_seqsum_basecall.log"
    benchmark:
        "{tmp}02_analysis_results/01_basecalling/{run}/nanocomp/MeBaPiNa_nanocomp_seqsum_basecall.benchmark.tsv"
    shell:
        "out1={output[0]}; "
        "cp ${{out1/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} ${{out1}}; "
        "out2={output[1]}; "
        "cp ${{out2/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} ${{out2}}; "
        "lg={log}; "
        "cat ${{lg/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} > {log}"

## NANOQC ##

rule plot_fastq_pipe_basecall:
    input:
        "{tmp}01_processed_data/01_basecalling/{run}/pass"
    output:
        temp("{tmp}02_analysis_results/01_basecalling/{run}/nanoqc/pipe.fastq") ## pipe didn't work (no fastq file extension) neighter did named pipes (started but nevenr finished)
    log:
        "{tmp}02_analysis_results/01_basecalling/{run}/nanoqc/MeBaPiNa_nanoqc_fastq_basecall.log"
    benchmark:
        "{tmp}02_analysis_results/01_basecalling/{run}/nanoqc/MeBaPiNa_nanoqc_fastq_basecall.benchmark.tsv"
    shell:
        "touch {output[0]}; "
        "lg={log}; "
        "cat ${{lg/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} > {log}"

rule plot_nanoqc_fastq_basecall:
    input:
        "{tmp}02_analysis_results/01_basecalling/{run}/nanoqc/pipe.fastq"
    output:
        "{tmp}02_analysis_results/01_basecalling/{run}/nanoqc/nanoQC.html"
    log:
        "{tmp}02_analysis_results/01_basecalling/{run}/nanoqc/MeBaPiNa_nanoqc_fastq_basecall.log"
    benchmark:
        "{tmp}02_analysis_results/01_basecalling/{run}/nanoqc/MeBaPiNa_nanoqc_fastq_basecall.benchmark.tsv"
    shell:
        "out1={output[0]}; "
        "cp ${{out1/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} ${{out1}}; "
        "lg={log}; "
        "cat ${{lg/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} > {log}"

## FASTQC ##

rule plot_fastqc_fastq_basecall:
    input:
        "{tmp}01_processed_data/01_basecalling/{run}/pass"
    output:
        "{tmp}02_analysis_results/01_basecalling/{run}/fastqc/stdin_fastqc.html"
    log:
        "{tmp}02_analysis_results/01_basecalling/{run}/fastqc/MeBaPiNa_fastqc_fastq_basecall.log"
    benchmark:
        "{tmp}02_analysis_results/01_basecalling/{run}/fastqc/MeBaPiNa_fastqc_fastq_basecall.benchmark.tsv"
    shell:
        "out1={output[0]}; "
        "cp ${{out1/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} ${{out1}}; "
        "lg={log}; "
        "cat ${{lg/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} > {log}"

## TRIM AND FILTER ##
#####################

## target rule for all output
def input_fastq(wildcards):
    from os import listdir
    ## get "pass" directory
    basecall_dir = checkpoints.basecall_raw.get(tmp=wildcards.tmp,run=wildcards.run).output[0]
    ## get barcode directory names within "pass" directory (excludes any barcodes without assigned reads)
    all_barc = listdir(basecall_dir)
    all_barc.sort()

    ## filtered fastq files for all barcodes
    input_list = ["{tmp}01_processed_data/02_trimming_filtering/{run}/" + barc + "/filtered.fastq" for barc in all_barc]
    ## return
    return input_list

## NANOPLOT ##

rule plot_nanoplot_fastq_filter:
    input:
        input_fastq
    output:
        "{tmp}02_analysis_results/02_trimming_filtering/{run}/nanoplot/NanoStats.txt",
        "{tmp}02_analysis_results/02_trimming_filtering/{run}/nanoplot/NanoPlot-report.html"
    log:
        "{tmp}02_analysis_results/02_trimming_filtering/{run}/nanoplot/MeBaPiNa_nanoplot_fastq_filter.log"
    benchmark:
        "{tmp}02_analysis_results/02_trimming_filtering/{run}/nanoplot/MeBaPiNa_nanoplot_fastq_filter.benchmark.tsv"
    shell:
        "out1={output[0]}; "
        "cp ${{out1/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} ${{out1}}; "
        "out2={output[1]}; "
        "cp ${{out2/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} ${{out2}}; "
        "lg={log}; "
        "cat ${{lg/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} > {log}"

## NANOCOMP ##

rule plot_nanocomp_fastq_filter:
    input:
        input_fastq
    output:
        "{tmp}02_analysis_results/02_trimming_filtering/{run}/nanocomp/NanoStats.txt",
        "{tmp}02_analysis_results/02_trimming_filtering/{run}/nanocomp/NanoComp-report.html"
    log:
        "{tmp}02_analysis_results/02_trimming_filtering/{run}/nanocomp/MeBaPiNa_nanocomp_fastq_filter.log"
    benchmark:
        "{tmp}02_analysis_results/02_trimming_filtering/{run}/nanocomp/MeBaPiNa_nanocomp_fastq_filter.benchmark.tsv"
    shell:
        "out1={output[0]}; "
        "cp ${{out1/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} ${{out1}}; "
        "out2={output[1]}; "
        "cp ${{out2/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} ${{out2}}; "
        "lg={log}; "
        "cat ${{lg/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} > {log}"

## NANOQC ##

rule plot_fastq_pipe_filter:
    input:
        input_fastq
    output:
        temp("{tmp}02_analysis_results/02_trimming_filtering/{run}/nanoqc/pipe.fastq") ## pipe (snakemake built-in or bash) didn't work neighter did named pipes (no fastq file extension)
    log:
        "{tmp}02_analysis_results/02_trimming_filtering/{run}/nanoqc/MeBaPiNa_fastq_pipe_filter.log"
    benchmark:
        "{tmp}02_analysis_results/02_trimming_filtering/{run}/nanoqc/MeBaPiNa_fastq_pipe_filter.benchmark.tsv"
    shell:
        "touch {output[0]}; "
        "lg={log}; "
        "cat ${{lg/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} > {log}"

rule plot_nanoqc_fastq_filter:
    input:
        "{tmp}02_analysis_results/02_trimming_filtering/{run}/nanoqc/pipe.fastq"
    output:
        "{tmp}02_analysis_results/02_trimming_filtering/{run}/nanoqc/nanoQC.html"
    log:
        "{tmp}02_analysis_results/02_trimming_filtering/{run}/nanoqc/MeBaPiNa_nanoqc_fastq_filter.log"
    benchmark:
        "{tmp}02_analysis_results/02_trimming_filtering/{run}/nanoqc/MeBaPiNa_nanoqc_fastq_filter.benchmark.tsv"
    shell:
        "out1={output[0]}; "
        "cp ${{out1/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} ${{out1}}; "
        "lg={log}; "
        "cat ${{lg/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} > {log}"

## FASTQC ##

rule plot_fastqc_fastq_filter:
    input:
        input_fastq
    output:
        "{tmp}02_analysis_results/02_trimming_filtering/{run}/fastqc/stdin_fastqc.html"
    log:
        "{tmp}02_analysis_results/02_trimming_filtering/{run}/fastqc/MeBaPiNa_fastqc_fastq_filter.log"
    benchmark:
        "{tmp}02_analysis_results/02_trimming_filtering/{run}/fastqc/MeBaPiNa_fastqc_fastq_filter.benchmark.tsv"
    shell:
        "out1={output[0]}; "
        "cp ${{out1/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} ${{out1}}; "
        "lg={log}; "
        "cat ${{lg/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} > {log}"

## OTU ##
#########

rule plot_qiime2_q2otupick:
    input:
        "{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/cluster_ftable.qza"
    output:
        "{tmp}02_analysis_results/03_otu_picking/{run}/{barc}/{reference}/q2otupick/index.html"
    log:
        "{tmp}02_analysis_results/03_otu_picking/{run}/{barc}/{reference}/q2otupick/MeBaPiNa_qiime2_q2otupick.log"
    benchmark:
        "{tmp}02_analysis_results/03_otu_picking/{run}/{barc}/{reference}/q2otupick/MeBaPiNa_qiime2_q2otupick.benchmark.tsv"
    shell:
        "out1={output[0]}; "
        "cp ${{out1/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} ${{out1}}; "
        "lg={log}; "
        "cat ${{lg/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} > {log}"

rule plot_qiime2_q2filter:
    input:
        "{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/filt_ftable.qza"
    output:
        "{tmp}02_analysis_results/03_otu_picking/{run}/{barc}/{reference}/q2filter/index.html"
    log:
        "{tmp}02_analysis_results/03_otu_picking/{run}/{barc}/{reference}/q2filter/MeBaPiNa_qiime2_q2filter.log"
    benchmark:
        "{tmp}02_analysis_results/03_otu_picking/{run}/{barc}/{reference}/q2filter/MeBaPiNa_qiime2_q2filter.benchmark.tsv"
    shell:
        "out1={output[0]}; "
        "cp ${{out1/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} ${{out1}}; "
        "lg={log}; "
        "cat ${{lg/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} > {log}"

rule plot_krona_q2rerep:
    input:
        output="{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}_{reftype}/filtered.kraken2",
        kronataxtab="{tmp}METADATA/Reference_Sequences/{reference}/krona/{reftype}/taxonomy.tab"
    output:
        "{tmp}02_analysis_results/03_otu_picking/{run}/{barc}/{reference}_{reftype}/krona.html"
    log:
        "{tmp}02_analysis_results/03_otu_picking/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_krona_q2rerep.log"
    benchmark:
        "{tmp}02_analysis_results/03_otu_picking/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_krona_q2rerep.benchmark.tsv"
    shell:
        "out1={output[0]}; "
        "cp ${{out1/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} ${{out1}}; "
        "lg={log}; "
        "cat ${{lg/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} > {log}"

## ALIGN ##
###########

rule plot_pycoqc_aligned:
    input:
        seqsum="{tmp}01_processed_data/01_basecalling/{run}/sequencing_summary/split", ## only folder is specified as output in splitting rule
        bam="{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/filteredsorted.bam"
    output:
        html="{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}/pycoqc.html",
        json="{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}/pycoqc.json"
    log:
        "{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}/MeBaPiNa_pycoqc_aligned.log"
    benchmark:
        "{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}/MeBaPiNa_pycoqc_aligned.benchmark.tsv"
    shell:
        "out1={output[0]}; "
        "cp ${{out1/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} ${{out1}}; "
        "out2={output[1]}; "
        "cp ${{out2/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} ${{out2}}; "
        "lg={log}; "
        "cat ${{lg/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} > {log}"

rule plot_refseq_coverage:
    input:
        "{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/refseq_coverage.tsv"
    output:
        covdist_plot="{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}/covdist.pdf",
        covpos_plot="{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}/covpos.pdf"
    shell:
        "out1={output[0]}; "
        "cp ${{out1/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} ${{out1}}; "
        "out2={output[1]}; "
        "cp ${{out2/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} ${{out2}}"

rule plot_krona_aligned_text:
    input:
        "{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}_{reftype}/aligned.counttaxlist"
    output:
        "{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}_{reftype}/krona.html"
    log:
        "{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_krona_aligned.log"
    benchmark:
        "{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_krona_aligned.benchmark.tsv"
    shell:
        "out1={output[0]}; "
        "cp ${{out1/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} ${{out1}}; "
        "lg={log}; "
        "cat ${{lg/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} > {log}"

## K-MER MAPPING ##
###################

rule plot_krona_kmermap_text:
    input:
        "{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/kmer.counttaxlist"
    output:
        "{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/krona_bracken.html"
    log:
        "{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_krona_kmermap_bracken.log"
    benchmark:
        "{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_krona_kmermap_bracken.benchmark.tsv"
    shell:
        "out1={output[0]}; "
        "cp ${{out1/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} ${{out1}}; "
        "lg={log}; "
        "cat ${{lg/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} > {log}"

rule plot_krona_kmermap_kraken:
    input:
        output="{tmp}01_processed_data/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/filtered.kraken2",
        kronataxtab="{tmp}METADATA/Reference_Sequences/{reference}/krona/{reftype}/taxonomy.tab"
    output:
        "{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/krona.html"
    log:
        "{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_krona_kmermap.log"
    benchmark:
        "{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_krona_kmermap.benchmark.tsv"
    shell:
        "out1={output[0]}; "
        "cp ${{out1/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} ${{out1}}; "
        "lg={log}; "
        "cat ${{lg/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} > {log}"
