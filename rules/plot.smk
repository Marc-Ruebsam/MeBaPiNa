##############
## PLOTTING ##
##############

## BASECALLING ##
#################

## NANOPLOT ##

rule plot_nanoplot_seqsum_basecall:
    input:
        "01_processed_data/01_basecalling/{run}/sequencing_summary.txt"
    output:
        "02_analysis_results/01_basecalling/{run}/nanoplot/NanoStats.txt"
    log:
        "02_analysis_results/01_basecalling/{run}/nanoplot/MeBaPiNa_nanoplot_seqsum_basecall.log"
    benchmark:
        "02_analysis_results/01_basecalling/{run}/nanoplot/MeBaPiNa_nanoplot_seqsum_basecall.benchmark.tsv"
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
        "--outdir 02_analysis_results/01_basecalling/{wildcards.run}/nanoplot "
        "--summary {input} > {log} 2>&1"

## PYCOQC ##

rule plot_pycoqc_seqsum_basecall:
    input:
        "01_processed_data/01_basecalling/{run}/sequencing_summary/sequencing_summary_downsampled.txt" #!# "01_processed_data/01_basecalling/{run}/sequencing_summary/sequencing_summary_sorted.txt"
    output:
        html="02_analysis_results/01_basecalling/{run}/pycoqc/pycoQC_report.html",
        json="02_analysis_results/01_basecalling/{run}/pycoqc/pycoQC_report.json"
    log:
        "02_analysis_results/01_basecalling/{run}/pycoqc/pycoqc_seqsum_basecall.log"
    benchmark:
        "02_analysis_results/01_basecalling/{run}/pycoqc/pycoqc_seqsum_basecall.benchmark.tsv"
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
        "01_processed_data/01_basecalling/{run}/sequencing_summary/sequencing_summary_sorted.txt"
    output:
        "02_analysis_results/01_basecalling/{run}/nanocomp/NanoStats.txt"
    log:
        "02_analysis_results/01_basecalling/{run}/nanocomp/MeBaPiNa_nanocomp_seqsum_basecall.log"
    benchmark:
        "02_analysis_results/01_basecalling/{run}/nanocomp/MeBaPiNa_nanocomp_seqsum_basecall.benchmark.tsv"
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
        "--outdir 02_analysis_results/01_basecalling/{wildcards.run}/nanocomp "
        "--summary {input} > {log} 2>&1"

## NANOQC ##

rule plot_fastq_pipe_basecall:
    input:
        "01_processed_data/01_basecalling/{run}/pass"
    output:
        temp("02_analysis_results/01_basecalling/{run}/nanoqc/pipe.fastq") ## pipe didn't work (no fastq file extension) neighter did named pipes (started but nevenr finished)
    log:
        "02_analysis_results/01_basecalling/{run}/nanoqc/MeBaPiNa_nanoqc_fastq_basecall.log"
    benchmark:
        "02_analysis_results/01_basecalling/{run}/nanoqc/MeBaPiNa_nanoqc_fastq_basecall.benchmark.tsv"
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
        "02_analysis_results/01_basecalling/{run}/nanoqc/pipe.fastq"
    output:
        "02_analysis_results/01_basecalling/{run}/nanoqc/nanoQC.html"
    log:
        "02_analysis_results/01_basecalling/{run}/nanoqc/MeBaPiNa_nanoqc_fastq_basecall.log"
    benchmark:
        "02_analysis_results/01_basecalling/{run}/nanoqc/MeBaPiNa_nanoqc_fastq_basecall.benchmark.tsv"
    conda:
        "../envs/nanopack.yml"
    params:
        "--minlen 240"
    shell:
        "nanoQC {params} "
        "--outdir 02_analysis_results/01_basecalling/{wildcards.run}/nanoqc "
        "{input} >> {log} 2>&1"

## FASTQC ##

rule plot_fastqc_fastq_basecall:
    input:
        "01_processed_data/01_basecalling/{run}/pass"
    output:
        "02_analysis_results/01_basecalling/{run}/fastqc/stdin_fastqc.html"
    log:
        "02_analysis_results/01_basecalling/{run}/fastqc/MeBaPiNa_fastqc_fastq_basecall.log"
    benchmark:
        "02_analysis_results/01_basecalling/{run}/fastqc/MeBaPiNa_fastqc_fastq_basecall.benchmark.tsv"
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
        "--outdir 02_analysis_results/01_basecalling/{wildcards.run}/fastqc "
        "stdin > {log} 2>&1"

## TRIM AND FILTER ##
#####################

## NANOPLOT ##

rule plot_nanoplot_fastq_filter:
    input:
        expand("01_processed_data/02_trimming_filtering/{run}/{barc}/filtered.fastq", run = RUNS, barc = SAMPLES.keys())
    output:
        "02_analysis_results/02_trimming_filtering/{run}/nanoplot/NanoStats.txt"
    log:
        "02_analysis_results/02_trimming_filtering/{run}/nanoplot/MeBaPiNa_nanoplot_fastq_filter.log"
    benchmark:
        "02_analysis_results/02_trimming_filtering/{run}/nanoplot/MeBaPiNa_nanoplot_fastq_filter.benchmark.tsv"
    conda:
        "../envs/nanopack.yml"
    threads:
        2
    params:
        "--plots kde hex dot",
        "--format svg",
        "--colormap viridis",
        "--color black", ## use NanoPlot --listcolors to get list of valid colors
        # "--raw", ## store "sequencing_summary"-like data
        "--downsample " + PLOT_SMPL, ## downlsampling
        "--verbose" ## or nothing to log
    shell:
        "NanoPlot --threads {threads} {params} "
        "--outdir 02_analysis_results/02_trimming_filtering/{wildcards.run}/nanoplot "
        "--fastq_rich {input} > {log} 2>&1"

## NANOCOMP ##

rule plot_nanocomp_fastq_filter:
    input:
        expand("01_processed_data/02_trimming_filtering/{run}/{barc}/filtered.fastq", run = RUNS, barc = SAMPLES.keys())
    output:
        "02_analysis_results/02_trimming_filtering/{run}/nanocomp/NanoStats.txt"
    log:
        "02_analysis_results/02_trimming_filtering/{run}/nanocomp/MeBaPiNa_nanocomp_fastq_filter.log"
    benchmark:
        "02_analysis_results/02_trimming_filtering/{run}/nanocomp/MeBaPiNa_nanocomp_fastq_filter.benchmark.tsv"
    conda:
        "../envs/nanopack.yml"
    threads:
        2 ## doesn't seem to use two cores
    params:
        "--maxlength " + PLOT_MAXLEN,
        "--plot violin", ## violin,box,ridge
        "--format svg",
        "--verbose", ## or nothing to log
        "--names " + " ".join(SAMPLES.values())
    shell:
        "NanoComp --threads {threads} {params} "
        "--outdir 02_analysis_results/02_trimming_filtering/{wildcards.run}/nanocomp "
        "--fastq {input} > {log} 2>&1"

## NANOQC ##

rule plot_fastq_pipe_filter:
    input:
        expand("01_processed_data/02_trimming_filtering/{run}/{barc}/filtered.fastq", run = RUNS, barc = SAMPLES.keys())
    output:
        temp("02_analysis_results/02_trimming_filtering/{run}/nanoqc/pipe.fastq") ## pipe (snakemake built-in or bash) didn't work neighter did named pipes (no fastq file extension)
    log:
        "02_analysis_results/02_trimming_filtering/{run}/nanoqc/MeBaPiNa_fastq_pipe_filter.log"
    benchmark:
        "02_analysis_results/02_trimming_filtering/{run}/nanoqc/MeBaPiNa_fastq_pipe_filter.benchmark.tsv"
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
        "02_analysis_results/02_trimming_filtering/{run}/nanoqc/pipe.fastq"
    output:
        "02_analysis_results/02_trimming_filtering/{run}/nanoqc/nanoQC.html"
    log:
        "02_analysis_results/02_trimming_filtering/{run}/nanoqc/MeBaPiNa_nanoqc_fastq_filter.log"
    benchmark:
        "02_analysis_results/02_trimming_filtering/{run}/nanoqc/MeBaPiNa_nanoqc_fastq_filter.benchmark.tsv"
    conda:
        "../envs/nanopack.yml"
    params:
        "--minlen 240"
    shell:
        "nanoQC {params} "
        "--outdir 02_analysis_results/02_trimming_filtering/{wildcards.run}/nanoqc "
        "{input} > {log} 2>&1"

## FASTQC ##

rule plot_fastqc_fastq_filter:
    input:
        expand("01_processed_data/02_trimming_filtering/{run}/{barc}/filtered.fastq", run = RUNS, barc = SAMPLES.keys())
    output:
        "02_analysis_results/02_trimming_filtering/{run}/fastqc/stdin_fastqc.html"
    log:
        "02_analysis_results/02_trimming_filtering/{run}/fastqc/MeBaPiNa_fastqc_fastq_filter.log"
    benchmark:
        "02_analysis_results/02_trimming_filtering/{run}/fastqc/MeBaPiNa_fastqc_fastq_filter.benchmark.tsv"
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
        "--outdir 02_analysis_results/02_trimming_filtering/{wildcards.run}/fastqc "
        "stdin > {log} 2>&1"

## ALIGN ##
###########

rule plot_pycoqc_bam_align:
    input:
        seqsum="01_processed_data/01_basecalling/{run}/sequencing_summary/split", ## only folder is specified as output in splitting rule
        bam="01_processed_data/03_alignment/{run}/{barc}/{altype}_sorted.bam"
    output:
        html="02_analysis_results/03_alignment/{run}/{barc}/pycoqc/{altype}.html",
        json="02_analysis_results/03_alignment/{run}/{barc}/pycoqc/{altype}.json"
    log:
        "02_analysis_results/03_alignment/{run}/pycoqc/MeBaPiNa_pycoqc_bam_align_{barc}_{altype}.log"
    benchmark:
        "02_analysis_results/03_alignment/{run}/pycoqc/MeBaPiNa_pycoqc_bam_align_{barc}_{altype}.benchmark.tsv"
    conda:
        "../envs/pycoqc.yml"
    params:
        ("--config Pineline/MeBaPiNa/scripts/pycoQC_config.json " if  ## use custom config (without coverage plot) for barcodes...
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

## K-MER MAPPING ##
###################

rule plot_krona_kraken2:
    input:
        output="01_processed_data/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/filtered.kraken2",
        reference="METADATA/Reference_Sequences/{reference}/krona/{reftype}"
    output:
        "02_analysis_results/03_kmer_mapping/{run}/{barc}/krona/{reference}_{reftype}/filtered.html"
    log:
        "02_analysis_results/03_kmer_mapping/{run}/{barc}/krona/{reference}_{reftype}/MeBaPiNa_krona_kraken2.log"
    benchmark:
                "02_analysis_results/03_kmer_mapping/{run}/{barc}/krona/{reference}_{reftype}/MeBaPiNa_krona_kraken2.benchmark.tsv"
    conda:
        "../envs/krona.yml"
    params:
        "-i", ## Include a wedge for queries with no hits.
        # "-q 2", ## Column of input files to use as query ID. Required if magnitude files are specified. [Default: '1']
        "-t 3" ## Column of input files to use as taxonomy ID. [Default: '2']
    shell:
        "ktImportTaxonomy {params} -tax {input.reference} {input.output} -o {output} > {log} 2>&1"

## CALIBRATION STRAIN ##
########################

rule plot_nanoplot_seqsum_calib:
    input:
        "01_processed_data/01_basecalling/{run}/sequencing_summary/split"
    output:
        "02_analysis_results/01_basecalling/{run}_calibration_strands/nanoplot/NanoStats.txt"
    log:
        "02_analysis_results/01_basecalling/{run}_calibration_strands/nanoplot/MeBaPiNa_nanoplot_seqsum_calib.log"
    benchmark:
        "02_analysis_results/01_basecalling/{run}_calibration_strands/nanoplot/MeBaPiNa_nanoplot_seqsum_calib.benchmark.tsv"
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
        "--outdir 02_analysis_results/01_basecalling/{wildcards.run}/nanoplot "
        "--summary {input}/lambda.txt > {log} 2>&1"

## alignment plots are created by the functions above
