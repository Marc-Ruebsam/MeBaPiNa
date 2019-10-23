rule nanoplot_seqsum:
    input:
        "01_basecalling/{run}/sequencing_summary.txt"
    output:
        "02_analysis/{run}/nanopack/nanoplot/LengthvsQualityScatterPlot_kde.svg"
    log:
        "02_analysis/{run}/nanopack/nanoplot/MeBaPiNa.log"
    benchmark:
        "02_analysis/{run}/nanopack/nanoplot/MeBaPiNa.benchmark.tsv"
    conda:
        "../envs/nanopack.yml"
    threads:
        config["machine"]["cpu"]
    params:
        ("--drop_outliers"),
        ("--plots kde hex dot"),
        ("--format svg"),
        ("--colormap plasma"),
        ("--color black"), ## use NanoPlot --listcolors to get list of valid colors
        ("--verbose"), ## or nothing to log
        ("--outdir 02_analysis/{run}/nanopack/nanoplot")
    shell:
        "NanoPlot --threads {threads} {params} --summary {input} 2> {log}"

rule nanoplot_fastq:
    input:
        "01_basecalling/{run}/pass"
    output:
        "02_analysis/{run}/nanopack/nanoplot/LengthvsQualityScatterPlot_kde.svg"
    log:
        "02_analysis/{run}/nanopack/nanoplot/MeBaPiNa.log"
    benchmark:
        "02_analysis/{run}/nanopack/nanoplot/MeBaPiNa.benchmark.tsv"
    conda:
        "../envs/nanopack.yml"
    threads:
        config["machine"]["cpu"]
    params:
        ("--drop_outliers"),
        ("--plots kde hex dot"),
        ("--format svg"),
        ("--colormap plasma"),
        ("--color black"), ## use NanoPlot --listcolors to get list of valid colors
        ("--verbose"), ## or nothing to log
        ("--outdir 02_analysis/{run}/nanopack/nanoplot")
    shell:
        "NanoPlot --threads {threads} {params} --fastq_rich {input}/* 2> {log}"

ruleorder: nanoplot_seqsum > nanoplot_fastq ## to solve disambiguities for now
