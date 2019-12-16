from snakemake.utils import min_version
## set minimum snakemake version
min_version("5.4")

## location of configuration and report specifications
configfile: "config.yaml"
report: "report/workflow.rst"

## set working directory to specified data location
workdir: config["experiment_directory"]["base"]

## define barcode and set default to "SQK-RAB204"
BAC_KIT = ("SQK-RAB204" if not config["guppy"]["bac_kit"] else config["guppy"]["bac_kit"]) ## sets default string if bac_kit is ""
## sample depth for downsampling in some plots
PLOT_SMPL = "100000"
## max read lengths in some plots
PLOT_MAXLEN = "4000"

## prevent unwanted extension of wildcards
wildcard_constraints:
    barc="\w+" ## is equalt to [a-zA-Z0-9_]+

## target output rule (the default end/output of the pipeline)
rule all:
    input:
        list(filter(None,[
        
        ## BASECALL ##
        
        ## general QC: all reads, including calibtation strads, intentional downsampling
        expand("02_analysis/{run}/basecall/nanoplot/NanoStats.txt",
        run=config["experiment_directory"]["run"]),
        ## general QC: all reads, forced downsampling
        expand("02_analysis/{run}/basecall/pycoqc/pycoQC_report.json",
        run=config["experiment_directory"]["run"]),
        ## per base QC: all reads, forced downsampling
        expand("02_analysis/{run}/basecall/nanoqc/nanoQC.html",
        run=config["experiment_directory"]["run"]),
        
        ## barcode QC: per barcode
        ("" if not config["guppy"]["bac_kit"] else ## "" if bac_kit is ""
        expand("02_analysis/{run}/basecall/nanocomp/NanoStats.txt",
        run=config["experiment_directory"]["run"])),
                
        ## calibration QC: only calinration strands
        ("" if not config["guppy"]["lam_DCS"] else ## "" if lam_DCS is False
        expand("02_analysis/{run}/basecall_calibration_strands/nanoplot/NanoStats.txt",
        run=config["experiment_directory"]["run"])),
        
        ## ALIGNMENT ##
        
        ## general QC: per barcode, intentional downsampling
        expand("02_analysis/{run}/align/MeBaPiNa_barcode_aggregation.txt",
        run=config["experiment_directory"]["run"]),
        
        ## calibration QC: only calinration strands
        ("" if not config["guppy"]["lam_DCS"] else ## "" if lam_DCS is False
        expand("02_analysis/{run}/align_calibration_strands/lambda_nanoplot/NanoStats.txt",
        run=config["experiment_directory"]["run"])),
        ("" if not config["guppy"]["lam_DCS"] else ## "" if lam_DCS is False
        expand("02_analysis/{run}/align_calibration_strands/lambda_pycoqc/pycoQC_report.json",
        run=config["experiment_directory"]["run"]))
        
        ]))

## include other rules
include: "rules/basecall.smk"
include: "rules/plot.smk"
include: "rules/align.smk"
include: "rules/misc.smk"
