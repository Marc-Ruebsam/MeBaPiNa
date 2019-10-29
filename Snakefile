from snakemake.utils import min_version
## set minimum snakemake version
min_version("5.2.0")

## location of configuration and report specifications
configfile: "config.yaml"
report: "report/workflow.rst"

## set working directory to specified data location
workdir: config["experiment_directory"]["base"]

## Allow users to fix the underlying OS via singularity.
# singularity: "docker://continuumio/miniconda3"

## target output rule (the default end/output of the pipeline)
rule all:
    input:
        list(filter(None,[
        ## evaluation of Lambda calibration strands only when specified
        (config["guppy"]["lam_DCS"] and
        expand("01_processeddata/{run}/align_lam/alignment.sam",
        run=config["experiment_directory"]["run"])),
        ## do NanoPlot of basecalled reads
        expand("02_analysis/{run}/basecall/nanoplot/NanoPlot-report.html",
        run=config["experiment_directory"]["run"]),
        ## align reads to reference
        expand("01_processeddata/{run}/align/alignment_sorted.bam", 
        run=config["experiment_directory"]["run"]),
        ## do NanoPlot of aligned reads reads
        expand("02_analysis/{run}/align/nanoplot/NanoPlot-report.html",
        run=config["experiment_directory"]["run"])
        ]))

## include other rules
include: "rules/basecall.smk"
include: "rules/plot.smk"
include: "rules/demultiplex.smk"
include: "rules/align.smk"
include: "rules/misc.smk"
