##############
## BASECALL ##
##############

rule nanoplot_seqsum:
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
        config["machine"]["cpu"]
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
        "--outdir 02_analysis/{wildcards.run}/basecall/nanoplot "
        "--summary {input} > {log} 2>&1"

rule pycoqc_seqsum:
    input:
        "01_processeddata/{run}/basecall/sequencing_summary.txt"
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
        "--sample " + PLOT_SMPL, ## downsampling
        "--verbose"
    shell:
        "pycoQC {params} "
        "--summary_file {input} "
        "--html_outfile {output.html} "
        "--json_outfile {output.json} > {log} 2>&1"

rule nanocomp_seqsum:
    input:
        "01_processeddata/{run}/basecall/sequencing_summary.txt"
    output:
        "02_analysis/{run}/basecall/nanocomp/NanoStats.txt"
    log:
        "02_analysis/{run}/basecall/nanocomp/MeBaPiNa_nanocomp_seqsum.log"
    benchmark:
        "02_analysis/{run}/basecall/nanocomp/MeBaPiNa_nanocomp_seqsum.benchmark.tsv"
    conda:
        "../envs/nanopack.yml"
    threads:
        config["machine"]["cpu"]
    params:
        "--maxlength 10000",
        "--barcoded",
        "--plot violin", ## violin,box,ridge
        "--format svg",
        "--verbose" ## or nothing to log
    shell:
        "NanoComp --threads {threads} {params} "
        "--outdir 02_analysis/{wildcards.run}/basecall/nanocomp "
        "--summary {input} > {log} 2>&1"

rule nanoqc:
    input:
        "02_analysis/{run}/basecall/nanoqc/pipe.fastq.gz"
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
        config["machine"]["cpu"]
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

###########
## ALIGN ##
###########

rule nanoplot_bam:
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
        config["machine"]["cpu"]
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

rule pycoqc_bam:
    input:
        seqsum="01_processeddata/{run}/basecall/sequencing_summary.txt",
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
        "--min_pass_qual 0",
        "--sample 100000",
        "--verbose"
    shell:
        "pycoQC {params} "
        "--summary_file {input.seqsum} "
        "--bam_file {input.bam} "
        "--html_outfile {output.html} "
        "--json_outfile {output.json} > {log} 2>&1"


# wub:
# "bam_accuracy.py -g bam_accuracy.tsv -l bam_accuracy_reads.tsv -r bam_accuracy.pdf -e ../../../01_processeddata/dummy_pool/align_calibration_strands/lambda_alignment_sorted.bam"
# "bam_alignment_length.py -t bam_alignment_length.tsv ../../../01_processeddata/dummy_pool/align_calibration_strands/lambda_alignment_sorted.bam"
# "bam_alignment_qc.py -f ../../../00_rawdata/reference_sequences/lambda/lambda_3.6kb.fasta -r bam_alignment_qc.pdf ../../../01_processeddata/dummy_pool/align_calibration_strands/lambda_alignment_sorted.bam"
#     "-x	Do not plot per-reference information. Default: False"
# "bam_count_reads.py -z ../../../00_rawdata/reference_sequences/lambda/lambda_3.6kb.fasta -g -t bam_count_reads.tsv ../../../01_processeddata/dummy_pool/align_calibration_strands/lambda_alignment_sorted.bam"
# "bam_gc_vs_qual.py -f ../../../00_rawdata/reference_sequences/lambda/lambda_3.6kb.fasta -r bam_gc_vs_qual.pdf -t bam_gc_vs_qual.tsv ../../../01_processeddata/dummy_pool/align_calibration_strands/lambda_alignment_sorted.bam"
# "bam_multi_qc -h"
# "bam_ref_base_coverage.py -f ../../../00_rawdata/reference_sequences/lambda/lambda_3.6kb.fasta -t bam_ref_base_coverage.tsv ../../../01_processeddata/dummy_pool/align_calibration_strands/lambda_alignment_sorted.bam"
# "bam_soft_clips_tab.py -t bam_soft_clips_tab.tsv ../../../01_processeddata/dummy_pool/align_calibration_strands/lambda_alignment_sorted.bam"
# "bias_explorer.py -r bias_explorer.pdf bam_count_reads.tsv"
