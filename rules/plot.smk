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
        "--drop_outliers",
        "--plots kde hex dot",
        "--format svg",
        "--colormap plasma",
        "--color black", ## use NanoPlot --listcolors to get list of valid colors
        "--verbose" ## or nothing to log
    shell:
        "NanoPlot --threads {threads} {params} "
        "--outdir 02_analysis/{wildcards.run}/basecall/nanoplot "
        "--summary {input} > {log} 2>&1"

rule pycoqc_seqsum:
    input:
        "01_processeddata/{run}/{align}/alignment_sorted.bam"
    output:
        html="02_analysis/{run}/{align}/pycoqc/pycoQC_report.html",
        json="02_analysis/{run}/{align}/pycoqc/pycoQC_report.json"
    log:
        "02_analysis/{run}/{align}/pycoqc/pycoQC_report.log"
    benchmark:
        "02_analysis/{run}/{align}/pycoqc/pycoQC_report.benchmark.tsv"
    conda:
        "../envs/pycoqc.yml"
    params:
        "--min_pass_qual 0",
        "--sample 100000",
        "--verbose"
    shell:
        "pycoQC {params} "
        "--summary_file {input} "
        "--html_outfile {output.html} "
        "--json_outfile {output.json} > {log} 2>&1"

# rule nanoplot_fastq:
#     input:
#         "01_processeddata/{run}/basecall/pass"
#         "01_processeddata/{run}/sequencing_summary.txt"
#     output:
#         "02_analysis/{run}/nanopack/nanoplot/LengthvsQualityScatterPlot_kde.svg"
#     log:
#         "02_analysis/{run}/nanopack/nanoplot/MeBaPiNa_nanoplot_fastq.log"
#     benchmark:
#         "02_analysis/{run}/nanopack/nanoplot/MeBaPiNa_nanoplot_fastq.benchmark.tsv"
#     conda:
#         "../envs/nanopack.yml"
#     threads:
#         config["machine"]["cpu"]
#     params:
#         "--drop_outliers",
#         "--plots kde hex dot",
#         "--format svg",
#         "--colormap plasma",
#         "--color black", ## use NanoPlot --listcolors to get list of valid colors
#         "--verbose" ## or nothing to log
#     shell:
#         "NanoPlot --threads {threads} {params} "
#         "--outdir 02_analysis/{run}/nanopack/nanoplot "
#         "--fastq_rich {input} > {log} 2>&1"

ruleorder: pycoqc_seqsum > nanoplot_seqsum # > nanoplot_fastq ## to solve disambiguities for now

rule nanoplot_bam:
    input:
        "01_processeddata/{run}/{align}/alignment_sorted.bam"
    output:
        "02_analysis/{run}/{align}/nanoplot/NanoStats.txt"
    log:
        "02_analysis/{run}/{align}/nanoplot/MeBaPiNa_nanoplot_bam.log"
    benchmark:
        "02_analysis/{run}/{align}/nanoplot/MeBaPiNa_nanoplot_bam.benchmark.tsv"
    conda:
        "../envs/nanopack.yml"
    threads:
        config["machine"]["cpu"]
    params:
        "--drop_outliers",
        "--plots kde hex dot",
        "--format svg",
        "--colormap plasma",
        "--color black", ## use NanoPlot --listcolors to get list of valid colors
        "--verbose" ## or nothing to log
    shell:
        "NanoPlot --threads {threads} {params} "
        "--outdir 02_analysis/{wildcards.run}/{wildcards.align}/nanoplot "
        "--bam {input} > {log} 2>&1"

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
    shell:
        "nanoQC "
        "--outdir 02_analysis/{wildcards.run}/basecall/nanoqc "
        "{input} > {log} 2>&1"
