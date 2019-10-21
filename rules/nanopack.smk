rule nanoplot:
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
        ("--plots kde hex"), #("--plots kde hex dot"),
        ("--format svg"),
        (lambda wildcards, threads: "--threads " + str(threads)), ## have to use function; function always has to have wildcards argument first
        ("--verbose"), ## or nothing to log
        ("--store"),
        ("--outdir 02_analysis/{run}/nanopack/nanoplot")
    shell:
        "NanoPlot {params} --summary {input} 2> {log}"

rule nanoplot_long:
    input:
        "01_basecalling/{run}/sequencing_summary.txt"
    output:
        "02_analysis/{run}/nanopack/nanoplot_long/LogTransformed_HistogramReadlength.svg"
    log:
        "02_analysis/{run}/nanopack/nanoplot_long/MeBaPiNa.log"
    benchmark:
        "02_analysis/{run}/nanopack/nanoplot_long/MeBaPiNa.benchmark.tsv"
    conda:
        "../envs/nanopack.yml"
    threads:
        config["machine"]["cpu"]
    params:
        ("--plots kde hex"), #("--plots kde hex dot"),
        ("--format svg"),
        (lambda wildcards, threads: "--threads " + str(threads)), ## have to use function; function always has to have wildcards argument first
        ("--verbose"), ## or nothing to log
        ("--store"),
        ("--outdir 02_analysis/{run}/nanopack/nanoplot_long")
    shell:
        "NanoPlot {params} --summary {input} 2> {log}"

rule nanoplot_bac:
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
        ("--plots kde hex"), #("--plots kde hex dot"),
        ("--format svg"),
        (lambda wildcards, threads: "--threads " + str(threads)), ## have to use function; function always has to have wildcards argument first
        ("--verbose"), ## or nothing to log
        ("--store"),
        ("--outdir 02_analysis/{run}/nanopack/nanoplot")
    shell:
        "NanoPlot {params} --summary {input} 2> {log}"

ruleorder: nanoplot > nanoplot_bac ## to solve disambiguities for now
