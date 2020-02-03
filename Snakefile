#####################################################
## MeBaPiNa - Meta Barcoding Pipeline for Nanopore ##
#####################################################

## author:      Marc Ruebsam (marc-ruebsam@hotmail.de)
## version:     0.7
## last change: 30.01.2020
## description: Pielpine for automated analysis of 16S metabarcoding samples.


## SNAKEMAKE CONFIGURATIONS ##
##############################

## import pakages
from snakemake.utils import min_version
import pandas as pd

## set minimum snakemake version
min_version("5.4")

## location of configuration and report specifications
configfile: "config.yaml"
report: "report/workflow.rst"

## set working directory to specified data location
workdir: config["experiments"]["project"]

## prevent unwanted extension of wildcards
wildcard_constraints:
    barc="[a-zA-Z0-9_-]+", 
    run="\w+", ## is equalt to [a-zA-Z0-9_]+
    reference="[a-zA-Z0-9_\.]+"

## METADATA AND THRESHOLDS ##
#############################

## load run information
METADATA = pd.read_excel( "METADATA/EXPERIMENT_SEQUENCING.xlsx", header = 1 )
## find meta data for required samples
METADATA = METADATA.loc[ METADATA['Sample name'].isin(config["experiments"]["samples"]), : ]
#!# currently only one single run analysis is supported. Use run with most sample overlaps
RUNS = METADATA['Run ID'].value_counts().sort_values(ascending=False).keys()[0]
METADATA = METADATA.loc[ METADATA['Run ID'].isin([RUNS]), : ]
## sample barcode information
# SAMPLES = METADATA['Sample name']
# BARCODES = METADATA['Barcode']
SAMPLES = pd.Series(METADATA['Sample name'].values,index=METADATA['Barcode']).to_dict()

## get run information
FLOWCELL = METADATA['Flow cell product'].unique()[0]
SEQ_KIT = METADATA['Sequencing kit'].unique()[0]
BAC_KIT = METADATA['Barcoding kit'].unique()[0]
LAM_DCS = METADATA['Lambda DCS'].unique()[0]
## set "SQK-RAB204" as default if bac_kit is ""
# BAC_KIT = ("SQK-RAB204" if not BAC_KIT else BAC_KIT)

## sample depth for downsampling in some plots
PLOT_SMPL = "100000"
## max read lengths in some plots
PLOT_MAXLEN = config["filtering"]["len_max"]


## PIELINE RULES AND END POINTS ##
##################################

## load rule set
include: "rules/basecall.smk"
include: "rules/plot.smk"
include: "rules/align.smk"
include: "rules/kmer_mapped.smk"
include: "rules/misc.smk"
include: "rules/report.smk"

## target output rule (the default end/output of the pipeline)
rule all:
    input:
        expand(list(filter(None,[

        ## BASECALL ##

        ## general QC: all reads, including calibtation strads, intentional downsampling
        "02_analysis_results/01_basecalling/{run}/nanoplot/NanoStats.txt",
        ## general QC: all reads, forced downsampling
        "02_analysis_results/01_basecalling/{run}/pycoqc/pycoQC_report.json",
        ## per base QC: all reads, forced downsampling
        "02_analysis_results/01_basecalling/{run}/nanoqc/nanoQC.html",
        ## read QC: all passed reads
        "02_analysis_results/01_basecalling/{run}/fastqc/stdin_fastqc.html",

        ## barcode QC: per barcode
        ("" if not BAC_KIT else ## "" if bac_kit is ""
        "02_analysis_results/01_basecalling/{run}/nanocomp/NanoStats.txt"),

        ## TRIM AND FILTER ##

        ## general QC: trimed and filtered barcoded reads, intentional downsampling
        "02_analysis_results/02_trimming_filtering/{run}/nanoplot/NanoStats.txt",
        ## per base QC: trimed and filtered barcoded reads, forced downsampling
        "02_analysis_results/02_trimming_filtering/{run}/nanoqc/nanoQC.html",
        ## read QC: trimed and filtered barcoded reads
        "02_analysis_results/02_trimming_filtering/{run}/fastqc/stdin_fastqc.html",

        ## barcode QC: trimed and filtered barcoded reads
        ("" if not BAC_KIT else ## "" if bac_kit is ""
        "02_analysis_results/02_trimming_filtering/{run}/nanocomp/NanoStats.txt"),

        ## ALIGNMENT ##

        ## general QC: per barcode, intentional downsampling
        "02_analysis_results/03_alignment/{run}/{barc}/pycoqc/filtered.html",
        # "02_analysis_results/03_kmer_mapping/{run}/bracken/barcode03/barcode03_S.bracken",

        ## CALIBRATION STRAIN ##

        ## calibration QC: only calinration strands
        ("" if not LAM_DCS else ## "" if lam_DCS is False
        "02_analysis_results/01_basecalling/{run}_calibration_strands/nanoplot/NanoStats.txt"),

        ## calibration QC: only calinration strands
        ("" if not LAM_DCS else ## "" if lam_DCS is False
        "02_analysis_results/03_alignment/{run}_calibration_strands/lambda_nanoplot/NanoStats.txt"),
        ("" if not LAM_DCS else ## "" if lam_DCS is False
        "02_analysis_results/03_alignment/{run}_calibration_strands/lambda_pycoqc/pycoQC_report.json")

        ])), run = RUNS, barc = SAMPLES.keys())
