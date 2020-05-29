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
    conda:
        "../envs/nanopack.yml"
    threads:
        2
    params:
        "--maxlength " + PLOT_MAXLEN, ## to keep it consistent with other plots
        "--drop_outliers",
        "--plots kde hex dot",
        "--format png",
        "--colormap viridis",
        "--color black", ## use NanoPlot --listcolors to get list of valid colors
        "--downsample " + PLOT_SMPL, ## downlsampling
        "--verbose" ## or nothing to log
    shell:
        "NanoPlot --threads {threads} {params} "
        "--outdir {wildcards.tmp}02_analysis_results/01_basecalling/{wildcards.run}/nanoplot "
        "--summary {input} > {log} 2>&1"

## PYCOQC ##

rule plot_pycoqc_seqsum_basecall:
    input:
        "{tmp}01_processed_data/01_basecalling/{run}/sequencing_summary/sequencing_summary_downsampled.txt" #!# "{tmp}01_processed_data/01_basecalling/{run}/sequencing_summary/sequencing_summary_sorted.txt"
    output:
        html="{tmp}02_analysis_results/01_basecalling/{run}/pycoqc/pycoQC_report.html",
        json="{tmp}02_analysis_results/01_basecalling/{run}/pycoqc/pycoQC_report.json"
    log:
        "{tmp}02_analysis_results/01_basecalling/{run}/pycoqc/pycoqc_seqsum_basecall.log"
    benchmark:
        "{tmp}02_analysis_results/01_basecalling/{run}/pycoqc/pycoqc_seqsum_basecall.benchmark.tsv"
    conda:
        "../envs/pycoqc.yml"
    params:
        "--min_pass_qual 0",
        "--filter_calibration", ## leave out calibration_strands
        #!# downsampling is done during generation of input "--sample " + PLOT_SMPL, ## downsampling
        "--min_barcode_percent 0.001", ## barcodes below 0.001% of reads are removed (1 in 100'000)
        "--verbose"
    shell:
        "pycoQC {params} "
        "--summary_file {input} "
        "--html_outfile {output.html} "
        "--json_outfile {output.json} > {log} 2>&1"

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
    conda:
        "../envs/nanopack.yml"
    threads:
        2 ## doesn't seem to use two cores
    params:
        "--maxlength " + PLOT_MAXLEN,
        "--barcoded",
        "--plot violin", ## violin,box,ridge
        "--format png",
        "--verbose" ## or nothing to log
    shell:
        "NanoComp --threads {threads} {params} "
        "--outdir {wildcards.tmp}02_analysis_results/01_basecalling/{wildcards.run}/nanocomp "
        "--summary {input} > {log} 2>&1"

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
        "find {input} -type f -name \"*.fastq\" -exec cat {{}} \\; | " ## concatenate all fastq files
        "awk -v seed=$RANDOM 'BEGIN{{ prnt=-4; nr=100; srand(seed) }}; " ## "falsify" print flag, set nr for fraction and set random seed. Note: downsampling is a fraction here not a number of reads
        "NR%4==1{{ " ## for every fourth line -> all headers
        "rhundr=1+int(rand()*nr); " ## get a random number between 1 and nr
        "if(rhundr==nr)" ## if the number is nr (by chance of 1/nr)
        "{{ prnt=NR }} " ## set prnt flag to current line number to...
        "}}; "
        "NR<prnt+4' " ## ...print this and the next three lines
        ">> {output} 2> {log}"

rule plot_nanoqc_fastq_basecall:
    input:
        "{tmp}02_analysis_results/01_basecalling/{run}/nanoqc/pipe.fastq"
    output:
        "{tmp}02_analysis_results/01_basecalling/{run}/nanoqc/nanoQC.html"
    log:
        "{tmp}02_analysis_results/01_basecalling/{run}/nanoqc/MeBaPiNa_nanoqc_fastq_basecall.log"
    benchmark:
        "{tmp}02_analysis_results/01_basecalling/{run}/nanoqc/MeBaPiNa_nanoqc_fastq_basecall.benchmark.tsv"
    conda:
        "../envs/nanopack.yml"
    params:
        "--minlen 240"
    shell:
        "nanoQC {params} "
        "--outdir {wildcards.tmp}02_analysis_results/01_basecalling/{wildcards.run}/nanoqc "
        "{input} >> {log} 2>&1"

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
    conda:
        "../envs/fastqc.yml"
    threads:
        2 ## might need adujustment if heap space is to low (java process stuck)
    params:
        "--min_length 240",
        "--format fastq"
    shell:
        "find {input} -type f -name \"*.fastq\" -exec cat {{}} \\; | "
        "fastqc --threads 8 {params} " #!!!#
        "--outdir {wildcards.tmp}02_analysis_results/01_basecalling/{wildcards.run}/fastqc "
        "stdin > {log} 2>&1"

## TRIM AND FILTER ##
#####################

## target rule for all output
def input_fastq(wildcards):
    from os import listdir
    ## get "pass" directory
    basecall_dir = checkpoints.basecall_raw.get(tmp=wildcards.tmp,run=wildcards.run).output[0]
    ## get barcode directory names within "pass" directory (excludes any barcodes without assigned reads)
    all_barc = listdir(basecall_dir)
    ## retain only barcodes containing one of the selected barcodes from the metadata (not unassigned)
    all_barc = [barc for barc in all_barc if barc in SAMPLES.keys()]

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
    conda:
        "../envs/nanopack.yml"
    threads:
        2
    params:
        "--plots kde hex dot",
        "--format png",
        "--colormap viridis",
        "--color black", ## use NanoPlot --listcolors to get list of valid colors
        # "--raw", ## store "sequencing_summary"-like data
        "--downsample " + PLOT_SMPL, ## downlsampling
        "--verbose" ## or nothing to log
    shell:
        "NanoPlot --threads {threads} {params} "
        "--outdir {wildcards.tmp}02_analysis_results/02_trimming_filtering/{wildcards.run}/nanoplot "
        "--fastq_rich {input} > {log} 2>&1"

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
    conda:
        "../envs/nanopack.yml"
    threads:
        2 ## doesn't seem to use two cores
    params:
        "--maxlength " + PLOT_MAXLEN,
        "--plot violin", ## violin,box,ridge
        "--format png",
        "--verbose", ## or nothing to log
        "--names " + " ".join(SAMPLES.keys()) ## join(SAMPLES.values()) works only with short names
    shell:
        "NanoComp --threads {threads} {params} "
        "--outdir {wildcards.tmp}02_analysis_results/02_trimming_filtering/{wildcards.run}/nanocomp "
        "--fastq {input} > {log} 2>&1"

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
        "cat {input} |" ## concatenate all fastq files
        "awk -v seed=$RANDOM 'BEGIN{{ prnt=-4; nr=100; srand(seed) }}; " ## "falsify" print flag, set nr for fraction and set random seed. Note: downsampling is a fraction here not a number of reads
        "NR%4==1{{ " #!# downsampling: for every fourth line -> all headers
        "rhundr=1+int(rand()*nr); " ## get a random number between 1 and nr
        "if(rhundr==nr)" ## if the number is nr (by chance of 1/nr)
        "{{ prnt=NR }} " ## set prnt flag to current line number to...
        "}}; "
        "NR<prnt+4' " ## ...print this and the next three lines
        ">> {output} 2> {log}"

rule plot_nanoqc_fastq_filter:
    input:
        "{tmp}02_analysis_results/02_trimming_filtering/{run}/nanoqc/pipe.fastq"
    output:
        "{tmp}02_analysis_results/02_trimming_filtering/{run}/nanoqc/nanoQC.html"
    log:
        "{tmp}02_analysis_results/02_trimming_filtering/{run}/nanoqc/MeBaPiNa_nanoqc_fastq_filter.log"
    benchmark:
        "{tmp}02_analysis_results/02_trimming_filtering/{run}/nanoqc/MeBaPiNa_nanoqc_fastq_filter.benchmark.tsv"
    conda:
        "../envs/nanopack.yml"
    params:
        "--minlen 240"
    shell:
        "nanoQC {params} "
        "--outdir {wildcards.tmp}02_analysis_results/02_trimming_filtering/{wildcards.run}/nanoqc "
        "{input} > {log} 2>&1"

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
    conda:
        "../envs/fastqc.yml"
    threads:
        2 ## might need adujustment if heap space is to low (java process stuck)
    params:
        "--min_length 240",
        "--format fastq"
    shell:
        "cat {input} | "
        "fastqc --threads 8 {params} " #!!!#
        "--outdir {wildcards.tmp}02_analysis_results/02_trimming_filtering/{wildcards.run}/fastqc "
        "stdin > {log} 2>&1"

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
    conda:
        "../envs/qiime2.yml"
    shell:
        "qiime feature-table summarize "
        "--i-table {input} "
        "--o-visualization {input}.qzv "
        "--verbose > {log} 2>&1; "
        "out_dir={output}; out_dir=\"${{out_dir/index.html/}}\" > {log} 2>&1; "
        "qiime tools export "
        "--input-path {input}.qzv "
        "--output-path ${{out_dir}} "
        ">> {log} 2>&1; "
        "rm {input}.qzv >> {log} 2>&1"

rule plot_qiime2_q2filter:
    input:
        "{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/filt_ftable.qza"
    output:
        "{tmp}02_analysis_results/03_otu_picking/{run}/{barc}/{reference}/q2filter/index.html"
    log:
        "{tmp}02_analysis_results/03_otu_picking/{run}/{barc}/{reference}/q2filter/MeBaPiNa_qiime2_q2filter.log"
    benchmark:
        "{tmp}02_analysis_results/03_otu_picking/{run}/{barc}/{reference}/q2filter/MeBaPiNa_qiime2_q2filter.benchmark.tsv"
    conda:
        "../envs/qiime2.yml"
    shell:
        "qiime feature-table summarize "
        "--i-table {input} "
        "--o-visualization {input}.qzv "
        "--verbose > {log} 2>&1; "
        "out_dir={output}; out_dir=\"${{out_dir/index.html/}}\" > {log} 2>&1; "
        "qiime tools export "
        "--input-path {input}.qzv "
        "--output-path ${{out_dir}} "
        ">> {log} 2>&1; "
        "rm {input}.qzv >> {log} 2>&1"

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
    conda:
        "../envs/krona.yml"
    params:
        "-i", ## Include a wedge for queries with no hits.
        # "-q 2", ## Column of input files to use as query ID. Required if magnitude files are specified. [Default: '1']
        "-t 3", ## Column of input files to use as taxonomy ID. [Default: '2']
        "-d 7" ## Maximum depth of wedges to include in the chart.
    shell:
        "reference={input.kronataxtab}; reference=\"${{reference/taxonomy.tab/}}\" > {log} 2>&1; "
        "ktImportTaxonomy {params} -tax  ${{reference}} {input.output} -o {output} > {log} 2>&1"

## for assignment with classifyer
# qiime metadata tabulate \
#   --m-input-file counttax.qza \
#   --o-visualization counttax.qzv

## ALIGN ##
###########

rule plot_pycoqc_aligned:
    input:
        seqsum="{tmp}01_processed_data/01_basecalling/{run}/sequencing_summary/split", ## only folder is specified as output in splitting rule
        bam="{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/filteredsorted.bam"
    output:
        html="{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}_{reftype}/pycoqc.html",
        json="{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}_{reftype}/pycoqc.json"
    log:
        "{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_pycoqc_aligned.log"
    benchmark:
        "{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_pycoqc_aligned.benchmark.tsv"
    conda:
        "../envs/pycoqc.yml"
    params:
        ("--config {tmp}Pipeline/MeBaPiNa/scripts/pycoQC_config.json " if  ## use custom config (without coverage plot) for barcodes...
        not config["reference"]["source"] == "zymobiomics" ## ...but not for the Zymo reference
        and not "{wildcards.barc}" == "lambda" else ""), ## or calibration strains
        "--min_pass_qual 0",
        "--sample " + PLOT_SMPL, ## downsampling
        "--verbose"
    shell:
        "pycoQC {params} "
        "--summary_file {input.seqsum}/{wildcards.barc}.txt " ## only folder is specified as output in splitting rule
        "--bam_file {input.bam} "
        "--html_outfile {output.html} "
        "--json_outfile {output.json} > {log} 2>&1"

rule plot_krona_aligned_text:
    input:
        "{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}_{reftype}/aligned.counttaxlist"
    output:
        "{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}_{reftype}/krona.html"
    log:
        "{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_krona_aligned.log"
    benchmark:
        "{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_krona_aligned.benchmark.tsv"
    conda:
        "../envs/krona.yml"
    params:
        "-n \"Root\"" ## name for highest taxon
    shell:
        "ktImportText {params} -o {output} {input} > {log} 2>&1"

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
    conda:
        "../envs/krona.yml"
    params:
        "-n \"Root\"" ## name for highest taxon
    shell:
        "ktImportText {params} -o {output} {input} > {log} 2>&1"

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
    conda:
        "../envs/krona.yml"
    params:
        "-i", ## Include a wedge for queries with no hits.
        # "-q 2", ## Column of input files to use as query ID. Required if magnitude files are specified. [Default: '1']
        "-t 3", ## Column of input files to use as taxonomy ID. [Default: '2']
        "-d 7" ## Maximum depth of wedges to include in the chart.
    shell:
        "reference={input.kronataxtab}; reference=\"${{reference/taxonomy.tab/}}\" > {log} 2>&1; "
        "ktImportTaxonomy {params} -tax  ${{reference}} {input.output} -o {output} > {log} 2>&1"

## CALIBRATION STRAIN ##
########################

# rule plot_nanoplot_seqsum_calib:
#     input:
#         "{tmp}01_processed_data/01_basecalling/{run}/sequencing_summary/split"
#     output:
#         "{tmp}02_analysis_results/01_basecalling/{run}_calibration_strands/nanoplot/NanoStats.txt"
#     log:
#         "{tmp}02_analysis_results/01_basecalling/{run}_calibration_strands/nanoplot/MeBaPiNa_nanoplot_seqsum_calib.log"
#     benchmark:
#         "{tmp}02_analysis_results/01_basecalling/{run}_calibration_strands/nanoplot/MeBaPiNa_nanoplot_seqsum_calib.benchmark.tsv"
#     conda:
#         "../envs/nanopack.yml"
#     threads:
#         2
#     params:
#         "--drop_outliers", ## other functions use "--maxlength 10000", to keep it consistent with other plots
#         "--plots kde hex dot",
#         "--format png",
#         "--colormap viridis",
#         "--color black", ## use NanoPlot --listcolors to get list of valid colors
#         "--downsample " + PLOT_SMPL, ## downlsampling
#         "--verbose" ## or nothing to log
#     shell:
#         "NanoPlot --threads {threads} {params} "
#         "--outdir {wildcards.tmp}02_analysis_results/01_basecalling/{wildcards.run}/nanoplot "
#         "--summary {input}/lambda.txt > {log} 2>&1"
#
# rule plot_pycoqc_bam_calib:
#     input:
#         seqsum="{tmp}01_processed_data/01_basecalling/{run}/sequencing_summary/split", ## only folder is specified as output in splitting rule
#         bam="{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/lambda/filteredsorted.bam"
#     output:
#         html="{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}_{reftype}/pycoqc/calibration.html",
#         json="{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}_{reftype}/pycoqc/calibration.json"
#     log:
#         "{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}_{reftype}/pycoqc/MeBaPiNa_pycoqc_bam_calib.log"
#     benchmark:
#         "{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}_{reftype}/pycoqc/MeBaPiNa_pycoqc_bam_calib.tsv"
#     conda:
#         "../envs/pycoqc.yml"
#     params:
#         ("--config {tmp}Pipeline/MeBaPiNa/scripts/pycoQC_config.json " if  ## use custom config (without coverage plot) for barcodes...
#         not config["reference"]["source"] == "zymobiomics" ## ...but not for the Zymo reference
#         and not "{wildcards.barc}" == "lambda" else ""), ## or calibration strains
#         "--min_pass_qual 0",
#         "--sample " + PLOT_SMPL, ## downsampling
#         "--verbose"
#     shell:
#         "pycoQC {params} "
#         "--summary_file {input.seqsum}/{wildcards.barc}.txt " ## only folder is specified as output in splitting rule
#         "--bam_file {input.bam} "
#         "--html_outfile {output.html} "
#         "--json_outfile {output.json} > {log} 2>&1"
