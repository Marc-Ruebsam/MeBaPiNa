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
METADATA = pd.read_excel( config["experiments"]["tmp"] + config["experiments"]["meta"], header = 1 )
## find meta data for required samples
METADATA = METADATA.loc[ METADATA['Sample name'].isin(config["experiments"]["samples"]), : ]
#!# currently only one single run analysis is supported. Use run with most sample overlaps
RUNS = [METADATA['Run ID'].value_counts().sort_values(ascending=False).keys()[0]]
METADATA = METADATA.loc[ METADATA['Run ID'].isin(RUNS), : ]
## sample barcode information
SAMPLES = pd.Series(METADATA['Sample name'].values,index=METADATA['Barcode']).to_dict()
TIMEPOINTS = pd.Series(METADATA['Zeitpunkt'].values,index=METADATA['Barcode']).to_dict()

## path to the csv file logs are written to
LOGS = config["experiments"]["tmp"] + config["experiments"]["meta"]

## get run information
FLOWCELL = METADATA['Flow cell product'].unique()[0]
SEQ_KIT = METADATA['Sequencing kit'].unique()[0]
BAC_KIT = METADATA['Barcoding kit'].unique()[0]
LAM_DCS = METADATA['Lambda DCS'].unique()[0]

## sample depth for downsampling in some plots
PLOT_SMPL = config["filtering"]["plot_sample"]
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

## target rule for all output
def input_all(wildcards):
    from os import listdir
    ## get "pass" directory
    basecall_dir = checkpoints.basecall_raw.get(tmp=config["experiments"]["tmp"],run=RUNS[0]).output[0]
    ## get barcode directory names within "pass" directory (excludes any barcodes without assigned reads)
    all_barc = listdir(basecall_dir)
    ## retain only barcodes containing one of the selected barcodes from the metadata (not unassigned)
    all_barc = [barc for barc in all_barc if barc in SAMPLES.keys()]

    ## report directories...
    all_barc_dir = (
    ## ...per PROMISE timepoint and sample as specified in the METADATA
    ["03_report/" + TPs + "/" + IDs + "/" + RUNs + "-" + barc + "/"
    for TPs,IDs,barc,RUNs in zip(TIMEPOINTS.values(), SAMPLES.values(), SAMPLES.keys(), METADATA['Run ID']) if ("PROM" in IDs) & (barc in all_barc)] +
    ## ...for all other sample specified in the METADATA
    ["03_report/" + "non-PROMISE_samples" + "/" + IDs + "/" + RUNs + "-" + barc + "/"
    for TPs,IDs,barc,RUNs in zip(TIMEPOINTS.values(), SAMPLES.values(), SAMPLES.keys(), METADATA['Run ID']) if (not "PROM" in IDs) & (barc in all_barc)])

    ## requested files...
    input_list = (
    ## REPORT ##
    ["{tmp}METADATA/ANALYSIS_PROGRESS_MANAGEMENT.csv"] +
    ## REFERENCE ##
    ["{tmp}03_report/Reference_Sequences/{reference}/reference_lengthdist.tsv", ## reference lenth distribution
    "{tmp}03_report/Reference_Sequences/{reference}/reference_lengthdist.pdf", ## reference lenth distribution
    "{tmp}03_report/Reference_Sequences/{reference}/reference_taxaranks.tsv"] + ## reference taxa distribution
    ## RAW READS ##
    ["{tmp}{barc_dir}read_base_counts.tsv"] + ## raw read statistics
    ## BASECALL ##
    ["{tmp}{barc_dir}01_basecalling-nanoplot-NanoPlot-report.html", ## general QC: all reads, including calibtation strads, intentional downsampling
    "{tmp}{barc_dir}01_basecalling-nanoplot-NanoStats.txt", ## general QC: all reads, including calibtation strads, intentional downsampling
    "{tmp}{barc_dir}01_basecalling-pycoqc-pycoQC_report.html", ## general QC: all reads, forced downsampling
    "{tmp}{barc_dir}01_basecalling-pycoqc-pycoQC_report.json", ## general QC: all reads, forced downsampling
    "{tmp}{barc_dir}01_basecalling-nanoqc-nanoQC.html", ## per base QC: all reads, forced downsampling
    "{tmp}{barc_dir}01_basecalling-fastqc-stdin_fastqc.html"] + ## read QC: all passed reads
    [x for x in
    ["{tmp}{barc_dir}01_basecalling-nanocomp-NanoComp-report.html", ## barcode QC: per barcode
    "{tmp}{barc_dir}01_basecalling-nanocomp-NanoStats.txt"] ## barcode QC: per barcode
    if BAC_KIT] + ## if BAC_KIT is not ""
    ## TRIM AND FILTER ##
    ["{tmp}{barc_dir}02_trimming_filtering-nanoplot-NanoPlot-report.html", ## general QC: trimed and filtered barcoded reads, intentional downsampling
    "{tmp}{barc_dir}02_trimming_filtering-nanoplot-NanoStats.txt", ## general QC: trimed and filtered barcoded reads, intentional downsampling
    "{tmp}{barc_dir}02_trimming_filtering-nanoqc-nanoQC.html", ## per base QC: trimed and filtered barcoded reads, forced downsampling
    "{tmp}{barc_dir}02_trimming_filtering-fastqc-stdin_fastqc.html"] + ## read QC: trimed and filtered barcoded reads
    [x for x in
    ["{tmp}{barc_dir}02_trimming_filtering-nanocomp-NanoComp-report.html", ## barcode QC: trimed and filtered barcoded reads
    "{tmp}{barc_dir}02_trimming_filtering-nanocomp-NanoStats.txt"] ## barcode QC: trimed and filtered barcoded reads
    if BAC_KIT] + ## if BAC_KIT is not ""
    ## OTU ##
    [x for x in
    ["{tmp}{barc_dir}03_otu_picking-{reference}-q2otupick-index.html", ## clustered reads
    "{tmp}{barc_dir}03_otu_picking-{reference}-q2filter-index.html", ## filtered reads
    "{tmp}{barc_dir}03_otu_picking-{reference}_{reftype}-krona.html", ## classified taxa
    "{tmp}{barc_dir}03_otu_picking-{reference}_{reftype}-kmer.counttaxlist", ## taxonomic classifications
    "{tmp}{barc_dir}otu_feature_counts-{reference}.tsv", ## feature statistics
    "{tmp}{barc_dir}otu_taxa_counts-{reference}_{reftype}.tsv"] ## taxa statistics
    if "otu" in config["methodologie"]] + ## if "otu" is selected
    ## ALIGNMENT ##
    [x for x in
    ["{tmp}{barc_dir}03_alignment-{reference}_{reftype}-pycoqc.html", ## per barcode, intentional downsampling
    "{tmp}{barc_dir}03_alignment-{reference}_{reftype}-pycoqc.json", ## per barcode, intentional downsampling
    "{tmp}{barc_dir}03_alignment-{reference}_{reftype}-krona.html", ## taxonomic classification
    "{tmp}{barc_dir}03_alignment-{reference}_{reftype}-aligned.counttaxlist", ## taxonomic classification
    "{tmp}{barc_dir}kmer_taxa_counts-{reference}_{reftype}.tsv", ## taxa statistics
    "{tmp}{barc_dir}kmer_retaxa_counts-{reference}_{reftype}.tsv"] ## taxa statistics after abundance reestimation
    if "align" in config["methodologie"]] + ## if "align" is selected
    ## K-MER MAPPING ##
    [x for x in
    ["{tmp}{barc_dir}03_kmer_mapping-{reference}_{reftype}-krona.html", ## taxonomic composition
    "{tmp}{barc_dir}03_kmer_mapping-{reference}_{reftype}-krona_bracken.html", ## taxonomic composition after reestimation
    "{tmp}{barc_dir}03_kmer_mapping-{reference}_{reftype}-kmer.counttaxlist"] ## taxonomic classification
    if "kmer" in config["methodologie"]] ) ## if "kmer" is selected

    ## expand for all barcodes
    input_list = expand(input_list,
    tmp = config["experiments"]["tmp"],
    barc_dir = all_barc_dir,
    reference = config['reference']['source'],
    reftype = config['reference']['rank'] )
    ## return
    return input_list

rule all_target:
    input:
        input_all
    shell:
        "awk 'NR == 1; NR > 1 {{print $0 | \"sort -n | uniq\"}}' \"{input[0]}\" > \"{input[0]}.temp\"; " ## store unique lines in temporary output
        "cat \"{input[0]}.temp\" > \"{input[0]}\"; rm \"{input[0]}.temp\"" ## convert temporary back to file
