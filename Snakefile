#######################################################################
## MeBaPiNa - Meta Barcoding Analysis Pipeline for Nanopore Datasets ##
#######################################################################

## author:      Marc Ruebsam (marc-ruebsam@hotmail.de)
## description: Pipeline for automated analysis of 16S metabarcoding samples.


## SNAKEMAKE CONFIGURATIONS ##
##############################

## import pakages
from snakemake.utils import min_version
import pandas as pd

## set minimum snakemake version
min_version("5.4")

## location of configuration and report specifications
configfile: "config.yaml"
report: "report/workflow.rst" #!# not yet implemented

## set working directory to specified data location
workdir: config["experiments"]["project"]

## prevent unwanted extension of wildcards
wildcard_constraints:
    barc="[a-zA-Z0-9]+",
    run="\w+", ## is equalt to [a-zA-Z0-9_]+
    reftype="[a-zA-Z0-9]+",
    reference="[a-zA-Z0-9]+"


## METADATA AND THRESHOLDS ##
#############################

## load run information
METADATA = pd.read_excel( config["experiments"]["meta"], header = 1 )
## find meta data for required samples
METADATA = METADATA.loc[ METADATA['Sample name'].isin(config["experiments"]["samples"]), : ]
#!# currently only one single run analysis is supported. Use run with most sample overlaps
RUNS = [METADATA['Run ID'].value_counts().sort_values(ascending=False).keys()[0]]
METADATA = METADATA.loc[ METADATA['Run ID'].isin(RUNS), : ]
## sample barcode information
SAMPLES = pd.Series(METADATA['Sample name'].values,index=METADATA['Barcode']).to_dict()
TIMEPOINTS = pd.Series(METADATA['Zeitpunkt'].values,index=METADATA['Barcode']).to_dict()

## get run information
FLOWCELL = METADATA['Flow cell product'].unique()[0]
SEQ_KIT = METADATA['Sequencing kit'].unique()[0]
BAC_KIT = METADATA['Barcoding kit'].unique()[0]
LAM_DCS = METADATA['Lambda DCS'].unique()[0]

## sample depth for downsampling in some plots
PLOT_SMPL = "100000"
## max read lengths in some plots
PLOT_MAXLEN = config["filtering"]["len_max"]


## PIELINE RULES ##
###################

## load rule set
include: "rules/basecall.smk"
include: "rules/align.smk"
include: "rules/kmer.smk"
include: "rules/otu.smk"
include: "rules/stats.smk"
include: "rules/plot.smk"
include: "rules/misc.smk"
include: "rules/report.smk"


## TARGET RULE ##
#################

## collection of all statistics (and few plots)
def input_stat(wildcards):
    from os import listdir
    ## get "pass" directory
    basecall_dir = checkpoints.basecall_raw.get(tmp=config["experiments"]["tmp"],run=RUNS[0]).output[0]
    ## get barcode directory names within "pass" directory
    all_barcs = listdir(basecall_dir)
    ## retain only folders containing one of the selected barcodes (not unassigned)
    all_barcs = [barc for barc in all_barcs if barc in SAMPLES.keys()]

    ## report directories per PROMISE timepoint and sample as specified in the METADATA
    promise_dirs = [config["experiments"]["tmp"] + "03_report/" + TPs + "/" + IDs + "/" + RUNs + "-" + barc + "/"
    for TPs,IDs,barc,RUNs in zip(TIMEPOINTS.values(), SAMPLES.values(), SAMPLES.keys(), METADATA['Run ID']) if ("PROM" in IDs) & (barc in all_barcs)]

    ## report directories for all other sample specified in the METADATA
    other_dirs = [config["experiments"]["tmp"] + "03_report/" + "non-PROMISE_samples" + "/" + IDs + "/" + RUNs + "-" + barc + "/"
    for TPs,IDs,barc,RUNs in zip(TIMEPOINTS.values(), SAMPLES.values(), SAMPLES.keys(), METADATA['Run ID']) if (not "PROM" in IDs) & (barc in all_barcs)]

    ## create file names with report dirs
    input_list = (
        ## RAW READS ##
        [stat_dir + "read_base_counts.tsv" for stat_dir in promise_dirs + other_dirs] +
        ## REFERENCE DATA ##
        ["{tmp}03_report/Reference_Sequences/{reference}/reference_lengthdist.tsv",
        "{tmp}03_report/Reference_Sequences/{reference}/reference_lengthdist.pdf",
        "{tmp}03_report/Reference_Sequences/{reference}/reference_taxaranks.tsv"] +
        ## REPORTS ##
        ["{tmp}METADATA/{run}-{reference}-{reftype}-reports.csv"]) ## from all_report rule
    input_list = expand(input_list,
    tmp = config["experiments"]["tmp"],
    run = RUNS,
    reference = config['reference']['source'],
    reftype = config['reference']['rank'] )
    ## return
    return input_list

## rule for stats
rule all:
    input:
        input_stat
