##############
## BASECALL ##
##############

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
        "--drop_outliers", ## other functions use "--maxlength 10000",
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
        "--maxlength " + PLOT_MAXLEN, # "--drop_outliers", ## keep consistent with other plots
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
        2
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
        "cat $(find {input} -type f -name \"*.fastq\") |" ## concatenate all fastq files
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

#####################
## TRIM AND FILTER ##
#####################

rule nanoplot_fastq_filter:
    input:
        "01_processeddata/{run}/filter/{barc}.fastq"
    output:
        "02_analysis/{run}/filter/{barc}_nanoplot/NanoStats.txt"
    log:
        "02_analysis/{run}/filter/{barc}_nanoplot/MeBaPiNa_nanoplot_seqsum.log"
    benchmark:
        "02_analysis/{run}/filter/{barc}_nanoplot/MeBaPiNa_nanoplot_seqsum.benchmark.tsv"
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
        "--outdir 02_analysis/{wildcards.run}/filter/{wildcards.barc}_nanoplot "
        "--fastq_rich {input} > {log} 2>&1"

rule fastq_pipe_filter:
    input:
        "01_processeddata/{run}/filter/{barc}.fastq"
    output:
        temp("02_analysis/{run}/filter/{barc}_nanoqc/pipe.fastq") ## pipe didn't work neighter did named pipes (no fastq file extension)
    log:
        "02_analysis/{run}/filter/{barc}_nanoqc/MeBaPiNa_fastq_pipe.log"
    benchmark:
        "02_analysis/{run}/filter/{barc}_nanoqc/MeBaPiNa_fastq_pipe.benchmark.tsv"
    shell:
        "awk -v seed=$RANDOM 'BEGIN{{ prnt=-4; nr=10; srand(seed) }}; " ## "falsify" print flag, set nr for fraction and set random seed. Note: downsampling is a fraction here not a number of reads
        "NR%4==1{{ " ## for every fourth line -> all headers
        "rhundr=1+int(rand()*nr); " ## get a random number between 1 and nr
        "if(rhundr==nr)" ## if the number is nr (by chance of 1/nr)
        "{{ prnt=NR }} " ## set prnt flag to current line number to... 
        "}}; "
        "NR<prnt+4' {input} " ## ...print this and the next three lines
        ">> {output} 2> {log}"

rule nanoqc_fastq_filter:
    input:
        "02_analysis/{run}/filter/{barc}_nanoqc/pipe.fastq"
    output:
        "02_analysis/{run}/filter/{barc}_nanoqc/nanoQC.html"
    log:
        "02_analysis/{run}/filter/{barc}_nanoqc/MeBaPiNa_nanoqc.log"
    benchmark:
        "02_analysis/{run}/filter/{barc}_nanoqc/MeBaPiNa_nanoqc.benchmark.tsv"
    conda:
        "../envs/nanopack.yml"
    params:
        "--minlen 240"
    shell:
        "nanoQC {params} "
        "--outdir 02_analysis/{wildcards.run}/filter/{wildcards.barc}_nanoqc "
        "{input} > {log} 2>&1"

###########
## ALIGN ##
###########

rule nanoplot_bam_align:
    input:
        bam="01_processeddata/{run}/{align}/{barc}_alignment_sorted.bam",
        bai="01_processeddata/{run}/{align}/{barc}_alignment_sorted.bam.bai"
    output:
        "02_analysis/{run}/{align}/{barc}_nanoplot/NanoStats.txt"
    log:
        "02_analysis/{run}/{align}/{barc}_nanoplot/MeBaPiNa_nanoplot_bam.log"
    benchmark:
        "02_analysis/{run}/{align}/{barc}_nanoplot/MeBaPiNa_nanoplot_bam.benchmark.tsv"
    conda:
        "../envs/nanopack.yml"
    threads:
        2
    params:
        "--maxlength 10000", ## to keep it consistent with other plots"--drop_outliers",
        "--alength", ## Use aligned read lengths rather than sequenced length (bam mode)
        "--plots kde hex dot",
        "--format svg",
        "--colormap plasma",
        "--color black", ## use NanoPlot --listcolors to get list of valid colors
        "--downsample " + PLOT_SMPL, ## downlsampling
        "--verbose" ## or nothing to log
    shell:
        "NanoPlot --threads {threads} {params} "
        "--outdir 02_analysis/{wildcards.run}/{wildcards.align}/{wildcards.barc}_nanoplot "
        "--bam {input.bam} > {log} 2>&1"

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
