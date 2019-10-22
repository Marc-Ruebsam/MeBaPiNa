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
        expand("01_basecalling/{run}/calibration_strands", 
        run=config["experiment_directory"]["run"])),
        ## demultiplexing only when specified
        (config["guppy"]["bac_kit"] and 
        expand("01_basecalling/{run}/pass_demultiplexed", 
        run=config["experiment_directory"]["run"])),
        ## do NanoPlot of basecalled reads
        expand("02_analysis/{run}/nanopack/nanoplot/LengthvsQualityScatterPlot_kde.svg", 
        run=config["experiment_directory"]["run"])   
        ]))

## include other rules        
include: "rules/other.smk"
include: "rules/guppy.smk"
include: "rules/nanopack.smk"
include: "rules/porechop.smk"
