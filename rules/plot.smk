rule nanoplot_seqsum:
    input:
        "01_processeddata/{run}/sequencing_summary.txt"
    output:
        "02_analysis/{run}/nanopack/nanoplot/LengthvsQualityScatterPlot_kde.svg"
    log:
        "02_analysis/{run}/nanopack/nanoplot/MeBaPiNa_nanoplot_seqsum.log"
    benchmark:
        "02_analysis/{run}/nanopack/nanoplot/MeBaPiNa_nanoplot_seqsum.benchmark.tsv"
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
        "--outdir 02_analysis/{run}/nanopack/nanoplot "
        "--summary {input} 2> {log}"

# rule nanoplot_fastq:
#     input:
#         "01_processeddata/{run}/pass"
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
#         "--fastq_rich {input} 2> {log}"

rule nanoqc:
    input:
        dummy_dependency="01_processeddata/{run}/sequencing_summary.txt",
        fastq="01_processeddata/{run}/pass"
    output:
        directory("02_analysis/{run}/nanopack/nanoqc")
    log:
        "02_analysis/{run}/nanopack/nanoqc/MeBaPiNa_nanoqc.log"
    benchmark:
        "02_analysis/{run}/nanopack/nanoqc/MeBaPiNa_nanoqc.benchmark.tsv"
    conda:
        "../envs/nanopack.yml"
    shell:
        "nanoQC --outdir 02_analysis/{run}/nanopack/nanoqc "
        "--fastq {input.fastq} 2> {log}"

# ruleorder: nanoplot_seqsum > nanoplot_fastq ## to solve disambiguities for now
