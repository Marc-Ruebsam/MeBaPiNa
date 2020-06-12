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
report: "report/workflow.rst" #!# not yet implemented

## set working directory to specified data location
workdir: config["experiments"]["project"]

## prevent unwanted extension of wildcards
wildcard_constraints:
    tmp="[a-zA-Z0-9_]+/",
    barc_dir="[a-zA-Z0-9_]+/",
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
## runs used to produce the samples
RUNS = [METADATA['Run ID'].value_counts().sort_values(ascending=False).keys()[0]] #!# only one run id can be used at a time
METADATA = METADATA.loc[ METADATA['Run ID'].isin(RUNS), : ]
## sample barcode information
SAMPLES = pd.Series(METADATA['Sample name'].values,index=METADATA['Barcode']).to_dict()
TIMEPOINTS = pd.Series(METADATA['Zeitpunkt'].values,index=METADATA['Barcode']).to_dict()

## get run information
FLOWCELL = METADATA['Flow cell product'].unique()[0]
SEQ_KIT = METADATA['Sequencing kit'].unique()[0]
BAC_KIT = METADATA['Barcoding kit'].unique()[0]
LAM_DCS = METADATA['Lambda DCS'].unique()[0]

## path to the csv file logs are written to
LOGS = config["experiments"]["tmp"] + config["experiments"]["log"]

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
    ## retain only barcodes containing one of the selected barcodes from the metadata (not unwanted barcodes)
    all_barc = [barc for barc in all_barc if barc in SAMPLES.keys()]
    all_barc.sort()

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
    ## LOGS ##
    ["{tmp}METADATA/ANALYSIS_PROGRESS_MANAGEMENT.csv"] +
    ## REFERENCE ##
    ["{tmp}03_report/Reference_Sequences/{reference}/reference_lengthdist.tsv", ## reference lenth distribution
    "{tmp}03_report/Reference_Sequences/{reference}/reference_lengthdist.pdf", ## reference lenth distribution
    "{tmp}03_report/Reference_Sequences/{reference}_{reftype}/reference_taxaranks.tsv"] + ## reference taxa distribution
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
    "{tmp}{barc_dir}03_otu_picking-{reference}_{reftype}-taxa_covdist.pdf", ## distribution of taxa abundance
    "{tmp}{barc_dir}03_otu_picking-{reference}-feature_counts.tsv", ## feature statistics
    "{tmp}{barc_dir}03_otu_picking-{reference}_{reftype}-taxa_counts.tsv", ## taxa statistics
    "{tmp}{barc_dir}03_otu_picking-{reference}_{reftype}-taxa_diversity.tsv"] ## diversity and richness measures
    if "otu" in config["methodologie"]] + ## if "otu" is selected
    ## ALIGNMENT ##
    [x for x in
    ["{tmp}{barc_dir}03_alignment-{reference}-pycoqc.html", ## per barcode, intentional downsampling
    "{tmp}{barc_dir}03_alignment-{reference}-pycoqc.json", ## per barcode, intentional downsampling
    "{tmp}{barc_dir}03_alignment-{reference}-covdist.pdf", ## per barcode, reads per reference sequence histogram
    "{tmp}{barc_dir}03_alignment-{reference}-covpos.pdf", ## per barcode, coverage over reference sequence positions
    "{tmp}{barc_dir}03_alignment-{reference}_{reftype}-krona.html", ## taxonomic classification
    "{tmp}{barc_dir}03_alignment-{reference}_{reftype}-aligned.counttaxlist", ## taxonomic classification
    "{tmp}{barc_dir}03_alignment-{reference}_{reftype}-taxa_covdist.pdf", ## distribution of taxa abundance
    "{tmp}{barc_dir}03_alignment-{reference}-alignment_rates.tsv", ## alignment statistics
    "{tmp}{barc_dir}03_alignment-{reference}_{reftype}-taxa_counts.tsv", ## taxa statistics
    "{tmp}{barc_dir}03_alignment-{reference}_{reftype}-taxa_diversity.tsv"] ## diversity and richness measures
    if "align" in config["methodologie"]] + ## if "align" is selected
    ## K-MER MAPPING ##
    [x for x in
    ["{tmp}{barc_dir}03_kmer_mapping-{reference}_{reftype}-krona.html", ## taxonomic composition
    "{tmp}{barc_dir}03_kmer_mapping-{reference}_{reftype}-krona_bracken.html", ## taxonomic composition after reestimation
    "{tmp}{barc_dir}03_kmer_mapping-{reference}_{reftype}-kmer.counttaxlist", ## taxonomic classification
    "{tmp}{barc_dir}03_kmer_mapping-{reference}_{reftype}-taxa_covdist.pdf", ## distribution of taxa abundance
    "{tmp}{barc_dir}03_kmer_mapping-{reference}_{reftype}-taxa_counts.tsv", ## taxa statistics
    "{tmp}{barc_dir}03_kmer_mapping-{reference}_{reftype}-retaxa_counts.tsv", ## taxa statistics after abundance reestimation
    "{tmp}{barc_dir}03_kmer_mapping-{reference}_{reftype}-taxa_diversity.tsv"] ## diversity and richness measures
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
        "echo \"Sample name;File/Directory;Action;Date;Checksum;Who;Description\" > \"{input[0]}.temp\"; " ## store header in temporary output
        "grep -v \"Sample name;File/Directory;Action;Date;Checksum;Who;Description\" \"{input[0]}\"  | sort -n | uniq >> \"{input[0]}.temp\"; " ## store unique lines in temporary output
        "cat \"{input[0]}.temp\" > \"{input[0]}\"; rm \"{input[0]}.temp\"" ## convert temporary back to file

## UPDATE REPORT ##
#################

## target rule for all output
def input_all(wildcards):
    from os import listdir
    ## get "pass" directory
    basecall_dir = checkpoints.basecall_raw.get(tmp=config["experiments"]["tmp"],run=RUNS[0]).output[0]
    ## get barcode directory names within "pass" directory (excludes any barcodes without assigned reads)
    all_barc = listdir(basecall_dir)
    ## retain only barcodes containing one of the selected barcodes from the metadata (not unwanted barcodes)
    all_barc = [barc for barc in all_barc if barc in SAMPLES.keys()]
    all_barc.sort()

    ## requested files...
    input_list = (
    ## LOGS ##
    ["{tmp}METADATA/ANALYSIS_PROGRESS_MANAGEMENT.csv"] +
    ## RAW READS ##
    ["{tmp}00_raw_data/{run}/MeBaPiNa_move_raw.report"] + ## REPORT
    ## BASECALL ##
    ["{tmp}00_raw_data/{run}/MeBaPiNa_basecall_raw_seqsum.report", ## REPORT
    "{tmp}00_raw_data/{run}/MeBaPiNa_basecall_raw_pass.report"] + ## REPORT
    ## TRIM AND FILTER ##
    ["{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/MeBaPiNa_trim_basecalled.report", ## REPORT
    "{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/MeBaPiNa_filter_trimmed.report"] + ## REPORT
    ## OTU ##
    [x for x in
    ["{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/MeBaPiNa_q2filter_uchime_ftable.report", ## REPORT
    "{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/MeBaPiNa_q2filter_uchime_centseq.report", ## REPORT
    "{tmp}02_analysis_results/03_otu_picking/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_kmermap_q2rereplicate.report", ## REPORT
    "{tmp}02_analysis_results/03_otu_picking/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_counttax_q2kmermap.report"] ## REPORT
    if "otu" in config["methodologie"]] + ## if "otu" is selected
    ## ALIGNMENT ##
    [x for x in
    ["{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/MeBaPiNa_filter_aligned.report",
    "{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_counttax_aligned.report"]
    if "align" in config["methodologie"]] + ## if "align" is selected
    ## K-MER MAPPING ##
    [x for x in
    ["{tmp}01_processed_data/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_kmermap_filtered.report", ## REPORT
    "{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_retax_kmermap.report", ## REPORT
    "{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_counttax_kmermap.report"] ## REPORT
    if "kmer" in config["methodologie"]] ) ## if "kmer" is selected

    ## expand for all barcodes
    input_list = expand(input_list,
    tmp = config["experiments"]["tmp"],
    run = RUNS,
    barc = all_barc,
    reference = config['reference']['source'],
    reftype = config['reference']['rank'] )
    ## return
    return input_list

rule update_report:
    input:
        input_all
    shell:
        "echo \"Sample name;File/Directory;Action;Date;Checksum;Who;Description\" > \"{input[0]}.temp\"; " ## store header in temporary output
        "grep -v \"Sample name;File/Directory;Action;Date;Checksum;Who;Description\" \"{input[0]}\"  | sort -n | uniq >> \"{input[0]}.temp\"; " ## store unique lines in temporary output
        "cat \"{input[0]}.temp\" > \"{input[0]}\"; rm \"{input[0]}.temp\"" ## convert temporary back to file
