##############
## BASECALL ##
##############

## NANOPLOT ##

rule nanoplot_seqsum_basecall:
    input:
        "01_processeddata/{run}/basecall/sequencing_summary.txt"
    output:
        "02_analysis/{run}/basecall/nanoplot/NanoStats.txt"
    log:
        "02_analysis/{run}/basecall/nanoplot/MeBaPiNa_nanoplot_seqsum.log"
    benchmark:
        "02_analysis/{run}/basecall/nanoplot/MeBaPiNa_nanoplot_seqsum.benchmark.tsv"
    conda:
        "../envs/nanopack.yml"
    threads:
        2
    params:
        "--maxlength " + PLOT_MAXLEN, ## to keep it consistent with other plots
        "--drop_outliers",
        "--plots kde hex dot",
        "--format svg",
        "--colormap viridis",
        "--color black", ## use NanoPlot --listcolors to get list of valid colors
        "--downsample " + PLOT_SMPL, ## downlsampling
        "--verbose" ## or nothing to log
    shell:
        "NanoPlot --threads {threads} {params} "
        "--outdir 02_analysis/{wildcards.run}/basecall/nanoplot "
        "--summary {input} > {log} 2>&1"

## PYCOQC ##

rule pycoqc_seqsum_basecall:
    input:
        "01_processeddata/{run}/basecall/sequencing_summary/sequencing_summary_downsampled.txt" # "01_processeddata/{run}/basecall/sequencing_summary/sequencing_summary_sorted.txt"
    output:
        html="02_analysis/{run}/basecall/pycoqc/pycoQC_report.html",
        json="02_analysis/{run}/basecall/pycoqc/pycoQC_report.json"
    log:
        "02_analysis/{run}/basecall/pycoqc/pycoQC_report.log"
    benchmark:
        "02_analysis/{run}/basecall/pycoqc/pycoQC_report.benchmark.tsv"
    conda:
        "../envs/pycoqc.yml"
    params:
        "--min_pass_qual 0",
        "--filter_calibration", ## leave out calibration_strands
        # downsampling is done during generation of input "--sample " + PLOT_SMPL, ## downsampling
        "--min_barcode_percent 0.001", ## barcodes below 0.001% of reads are removed (1 in 100'000)
        "--verbose"
    shell:
        "pycoQC {params} "
        "--summary_file {input} "
        "--html_outfile {output.html} "
        "--json_outfile {output.json} > {log} 2>&1"

## NANOCOMP ##

rule nanocomp_seqsum_basecall:
    input:
        "01_processeddata/{run}/basecall/sequencing_summary/sequencing_summary_sorted.txt"
    output:
        "02_analysis/{run}/basecall/nanocomp/NanoStats.txt"
    log:
        "02_analysis/{run}/basecall/nanocomp/MeBaPiNa_nanocomp_seqsum.log"
    benchmark:
        "02_analysis/{run}/basecall/nanocomp/MeBaPiNa_nanocomp_seqsum.benchmark.tsv"
    conda:
        "../envs/nanopack.yml"
    threads:
        2 ## doesn't seem to use two cores
    params:
        "--maxlength " + PLOT_MAXLEN,
        "--barcoded",
        "--plot violin", ## violin,box,ridge
        "--format svg",
        "--verbose" ## or nothing to log
    shell:
        "NanoComp --threads {threads} {params} "
        "--outdir 02_analysis/{wildcards.run}/basecall/nanocomp "
        "--summary {input} > {log} 2>&1"

## NANOQC ##

rule fastq_pipe_basecall:
    input:
        "01_processeddata/{run}/basecall/pass"
    output:
        temp("02_analysis/{run}/basecall/nanoqc/pipe.fastq") ## pipe didn't work neighter did named pipes (no fastq file extension)
    log:
        "02_analysis/{run}/basecall/nanoqc/MeBaPiNa_fastq_pipe.log"
    benchmark:
        "02_analysis/{run}/basecall/nanoqc/MeBaPiNa_fastq_pipe.benchmark.tsv"
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

rule nanoqc_fastq_basecall:
    input:
        "02_analysis/{run}/basecall/nanoqc/pipe.fastq"
    output:
        "02_analysis/{run}/basecall/nanoqc/nanoQC.html"
    log:
        "02_analysis/{run}/basecall/nanoqc/MeBaPiNa_nanoqc.log"
    benchmark:
        "02_analysis/{run}/basecall/nanoqc/MeBaPiNa_nanoqc.benchmark.tsv"
    conda:
        "../envs/nanopack.yml"
    params:
        "--minlen 240"
    shell:
        "nanoQC {params} "
        "--outdir 02_analysis/{wildcards.run}/basecall/nanoqc "
        "{input} > {log} 2>&1"

## FASTQC ##

rule fastqc_fastq_basecall:
    input:
        "01_processeddata/{run}/basecall/pass"
    output:
        "02_analysis/{run}/basecall/fastqc/stdin_fastqc.html"
    log:
        "02_analysis/{run}/basecall/fastqc/MeBaPiNa_fastqc.log"
    benchmark:
        "02_analysis/{run}/basecall/fastqc/MeBaPiNa_fastqc.benchmark.tsv"
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
        "--outdir 02_analysis/{wildcards.run}/basecall/fastqc "
        "stdin > {log} 2>&1"

#####################
## TRIM AND FILTER ##
#####################

## INPUT FUNCTION ##

def input_fastq_filter(wildcards):
    from os import listdir
    ## "pass" directory
    basecall_dir = checkpoints.guppy.get(run=wildcards.run).output[0]
    ## directory names within "pass" directory
    barcs = listdir(basecall_dir)
    ## retain only strings containing "barcode"
    barcs = [i for i in barcs if "barcode" in i]
    ## create file names with barcodes
    barc_input = expand("01_processeddata/{run}/filter/{barc}.fastq",
        run=wildcards.run,
        barc=barcs)
    barc_input.sort()
    return barc_input

## NANOPLOT ##

rule nanoplot_fastq_filter:
    input:
        input_fastq_filter
    output:
        "02_analysis/{run}/filter/nanoplot/NanoStats.txt"
    log:
        "02_analysis/{run}/filter/nanoplot/MeBaPiNa_nanoplot_seqsum.log"
    benchmark:
        "02_analysis/{run}/filter/nanoplot/MeBaPiNa_nanoplot_seqsum.benchmark.tsv"
    conda:
        "../envs/nanopack.yml"
    threads:
        2
    params:
        "--plots kde hex dot",
        "--format svg",
        "--colormap viridis",
        "--color black", ## use NanoPlot --listcolors to get list of valid colors
        "--downsample " + PLOT_SMPL, ## downlsampling
        "--verbose" ## or nothing to log
    shell:
        "NanoPlot --threads {threads} {params} "
        "--outdir 02_analysis/{wildcards.run}/filter/nanoplot "
        "--fastq_rich {input} > {log} 2>&1"

## NANOCOMP ##

rule nanocomp_fastq_filter:
    input:
        input_fastq_filter
    output:
        "02_analysis/{run}/filter/nanocomp/NanoStats.txt"
    log:
        "02_analysis/{run}/filter/nanocomp/MeBaPiNa_nanocomp_seqsum.log"
    benchmark:
        "02_analysis/{run}/filter/nanocomp/MeBaPiNa_nanocomp_seqsum.benchmark.tsv"
    conda:
        "../envs/nanopack.yml"
    threads:
        2 ## doesn't seem to use two cores
    params:
        "--maxlength " + PLOT_MAXLEN,
        "--plot violin", ## violin,box,ridge
        "--format svg",
        "--verbose" ## or nothing to log
    shell:
        "NanoComp --threads {threads} {params} "
        "--outdir 02_analysis/{wildcards.run}/filter/nanocomp "
        "--names $(echo {input} | sed 's#01_processeddata/{wildcards.run}/filter/##g' | sed 's#.fastq##g') "
        "--fastq {input} > {log} 2>&1"

## NANOQC ##

rule fastq_pipe_filter:
    input:
        input_fastq_filter
    output:
        temp("02_analysis/{run}/filter/nanoqc/pipe.fastq") ## pipe didn't work neighter did named pipes (no fastq file extension)
    log:
        "02_analysis/{run}/filter/nanoqc/MeBaPiNa_fastq_pipe.log"
    benchmark:
        "02_analysis/{run}/filter/nanoqc/MeBaPiNa_fastq_pipe.benchmark.tsv"
    shell:
        "cat {input} |" ## concatenate all fastq files
        "awk -v seed=$RANDOM 'BEGIN{{ prnt=-4; nr=100; srand(seed) }}; " ## "falsify" print flag, set nr for fraction and set random seed. Note: downsampling is a fraction here not a number of reads
        "NR%4==1{{ " ## for every fourth line -> all headers
        "rhundr=1+int(rand()*nr); " ## get a random number between 1 and nr
        "if(rhundr==nr)" ## if the number is nr (by chance of 1/nr)
        "{{ prnt=NR }} " ## set prnt flag to current line number to... 
        "}}; "
        "NR<prnt+4' " ## ...print this and the next three lines
        ">> {output} 2> {log}"

rule nanoqc_fastq_filter:
    input:
        "02_analysis/{run}/filter/nanoqc/pipe.fastq"
    output:
        "02_analysis/{run}/filter/nanoqc/nanoQC.html"
    log:
        "02_analysis/{run}/filter/nanoqc/MeBaPiNa_nanoqc.log"
    benchmark:
        "02_analysis/{run}/filter/nanoqc/MeBaPiNa_nanoqc.benchmark.tsv"
    conda:
        "../envs/nanopack.yml"
    params:
        "--minlen 240"
    shell:
        "nanoQC {params} "
        "--outdir 02_analysis/{wildcards.run}/filter/nanoqc "
        "{input} > {log} 2>&1"

## FASTQC ##

rule fastqc_fastq_filter:
    input:
        input_fastq_filter
    output:
        "02_analysis/{run}/filter/fastqc/stdin_fastqc.html"
    log:
        "02_analysis/{run}/filter/fastqc/MeBaPiNa_fastqc.log"
    benchmark:
        "02_analysis/{run}/filter/fastqc/MeBaPiNa_fastqc.benchmark.tsv"
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
        "--outdir 02_analysis/{wildcards.run}/filter/fastqc "
        "stdin > {log} 2>&1"

###########
## ALIGN ##
###########

# rule nanoplot_bam_align:
#     input:
#         bam="01_processeddata/{run}/{align}/{barc}_alignment_sorted.bam",
#         bai="01_processeddata/{run}/{align}/{barc}_alignment_sorted.bam.bai"
#     output:
#         "02_analysis/{run}/{align}/{barc}_nanoplot/NanoStats.txt"
#     log:
#         "02_analysis/{run}/{align}/{barc}_nanoplot/MeBaPiNa_nanoplot_bam.log"
#     benchmark:
#         "02_analysis/{run}/{align}/{barc}_nanoplot/MeBaPiNa_nanoplot_bam.benchmark.tsv"
#     conda:
#         "../envs/nanopack.yml"
#     threads:
#         2
#     params:
#         "--maxlength " + PLOT_MAXLEN, ## to keep it consistent with other plots
#         "--drop_outliers",
#         "--alength", ## Use aligned read lengths rather than sequenced length (bam mode)
#         "--plots kde hex dot",
#         "--format svg",
#         "--colormap plasma",
#         "--color black", ## use NanoPlot --listcolors to get list of valid colors
#         "--downsample " + PLOT_SMPL, ## downlsampling
#         "--verbose" ## or nothing to log
#     shell:
#         "NanoPlot --threads {threads} {params} "
#         "--outdir 02_analysis/{wildcards.run}/{wildcards.align}/{wildcards.barc}_nanoplot "
#         "--bam {input.bam} > {log} 2>&1"

rule pycoqc_bam_align:
    input:
        seqsum="01_processeddata/{run}/basecall/sequencing_summary/barcode", ## only folder is specified as output in splitting rule
        bam="01_processeddata/{run}/{align}/{barc}_alignment_sorted.bam"
    output:
        html="02_analysis/{run}/{align}/{barc}_pycoqc/pycoQC_report.html",
        json="02_analysis/{run}/{align}/{barc}_pycoqc/pycoQC_report.json"
    log:
        "02_analysis/{run}/{align}/{barc}_pycoqc/pycoQC_report.log"
    benchmark:
        "02_analysis/{run}/{align}/{barc}_pycoqc/pycoQC_report.benchmark.tsv"
    conda:
        "../envs/pycoqc.yml"
    params:
        ("--config MeBaPiNa/scripts/pycoQC_config.json " if  ## use custom config (without coverage plot) for barcodes...
        not config["align"]["reference"] == "ZymoBIOMICS.STD.refseq.v2/ssrRNAs/_all" ## ...but not for the Zymo reference
        and not "{wildcards.barc}" == "lambda" else ""), ## or calibration strains
        "--min_pass_qual 0",
        "--sample " + PLOT_SMPL, ## downsampling
        "--verbose"
    shell:
        "pycoQC {params} "
        "--summary_file {input.seqsum}/sequencing_summary_{wildcards.barc}.txt " ## only folder is specified as output in splitting rule
        "--bam_file {input.bam} "
        "--html_outfile {output.html} "
        "--json_outfile {output.json} > {log} 2>&1"

########################
## CALIBRATION STRAIN ##
########################

rule nanoplot_fastq_calib:
    input:
        "01_processeddata/{run}/basecall/calibration_strands"
    output:
        "02_analysis/{run}/basecall_calibration_strands/nanoplot/NanoStats.txt"
    log:
        "02_analysis/{run}/basecall_calibration_strands/nanoplot/MeBaPiNa_nanoplot_seqsum.log"
    benchmark:
        "02_analysis/{run}/basecall_calibration_strands/nanoplot/MeBaPiNa_nanoplot_seqsum.benchmark.tsv"
    conda:
        "../envs/nanopack.yml"
    threads:
        2
    params:
        "--drop_outliers", ## other functions use "--maxlength 10000", to keep it consistent with other plots
        "--plots kde hex dot",
        "--format svg",
        "--colormap viridis",
        "--color black", ## use NanoPlot --listcolors to get list of valid colors
        "--downsample " + PLOT_SMPL, ## downlsampling
        "--verbose" ## or nothing to log
    shell:
        "NanoPlot --threads {threads} {params} "
        "--outdir 02_analysis/{wildcards.run}/basecall_calibration_strands/nanoplot "
        "--fastq_rich {input}/* > {log} 2>&1"

## alignment plots are created by the functions above
