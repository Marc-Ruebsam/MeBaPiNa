from snakemake.utils import min_version
## set minimum snakemake version
min_version("5.4")

## location of configuration and report specifications
configfile: "config.yaml"
report: "report/workflow.rst"

## set working directory to specified data location
workdir: config["experiment_directory"]["base"]

## define barcode and set default to "SQK-RAB204"
BAC_KIT = ("SQK-RAB204" if config["guppy"]["bac_kit"] == "" else config["guppy"]["bac_kit"])
## sample depth for plotiing
PLOT_SMPL = "100000"

## Allow users to fix the underlying OS via singularity.
# singularity: "docker://continuumio/miniconda3"

## prevent unwanted extension of wildcards
wildcard_constraints:
    barc="\w+" ## is equalt to [a-zA-Z0-9_]+

## target output rule (the default end/output of the pipeline)
rule all:
    input:
        list(filter(None,[
        ## BASECALL ##
        ## comprehensive QC of all reads INcluding calibtation strands
        expand("02_analysis/{run}/basecall/nanoplot/NanoStats.txt",
        run=config["experiment_directory"]["run"]),
        ## comprehensive QC of all reads EXcluding calibrations strands
        expand("02_analysis/{run}/basecall/pycoqc/pycoQC_report.json",
        run=config["experiment_directory"]["run"]),
        ## comparison of the barcodes INcluding calibtation strands
        expand("02_analysis/{run}/basecall/nanocomp/NanoStats.txt",
        run=config["experiment_directory"]["run"]),
        ## per base quality overview EXcluding calibtation strands
        expand("02_analysis/{run}/basecall/nanoqc/nanoQC.html",
        run=config["experiment_directory"]["run"]),
        ## calibration strand specific comprehensive QC
        (config["guppy"]["lam_DCS"] and
        expand("02_analysis/{run}/basecall_calibration_strands/nanoplot/NanoStats.txt",
        run=config["experiment_directory"]["run"]))
        ]))

## include other rules
include: "rules/basecall.smk"
include: "rules/plot.smk"
include: "rules/demultiplex.smk"
include: "rules/align.smk"
include: "rules/misc.smk"
