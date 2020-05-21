#####################################################
## MeBaPiNa - Meta Barcoding Pipeline for Nanopore ##
#####################################################

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
report: "report/workflow.rst"

## set working directory to specified data location
workdir: config["experiments"]["project"]

## prevent unwanted extension of wildcards
wildcard_constraints:
    barc="[a-zA-Z0-9_-]+",
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
RUNS = METADATA['Run ID'].value_counts().sort_values(ascending=False).keys()[0]
METADATA = METADATA.loc[ METADATA['Run ID'].isin([RUNS]), : ]
## sample barcode information
SAMPLES = pd.Series(METADATA['Sample name'].values,index=METADATA['Barcode']).to_dict()

## get run information
FLOWCELL = METADATA['Flow cell product'].unique()[0]
SEQ_KIT = METADATA['Sequencing kit'].unique()[0]
BAC_KIT = METADATA['Barcoding kit'].unique()[0]
LAM_DCS = METADATA['Lambda DCS'].unique()[0]

## sample depth for downsampling in some plots
PLOT_SMPL = "100000"
## max read lengths in some plots
PLOT_MAXLEN = config["filtering"]["len_max"]


## PIELINE RULES AND END POINTS ##
##################################

## load rule set
include: "rules/basecall.smk"
include: "rules/align.smk"
include: "rules/kmer.smk"
include: "rules/otu.smk"
include: "rules/stats.smk"
include: "rules/plot.smk"
include: "rules/misc.smk"
include: "rules/report.smk"

## target output rule (the default end/output of the pipeline)
def input_barc(wildcards):
    from os import listdir
    ## get "pass" directory
    basecall_dir = checkpoints.basecalling_raw.get(tmp=config["experiments"]["tmp"],run=RUNS).output[1]
    ## get barcode directory names within "pass" directory
    all_barcs = listdir(basecall_dir)
    ## retain only folders containing one of the selected barcodes
    all_barcs = [barc for barc in all_barcs if barc in SAMPLES.keys()]
    ## create file names with barcodes
    input_list = expand(list(filter(None,

        ## BASECALL ##

        ["{tmp}00_raw_data/{run}/MeBaPiNa_moving_raw.report",
        "{tmp}00_raw_data/{run}/MeBaPiNa_basecalling_raw.report"] +

        ## TRIM AND FILTER ##

        ["{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/MeBaPiNa_trimming_basecalled.report",
        "{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/MeBaPiNa_filtering_trimmed.report"] +

        ## OTU ##

        ([""] if not "otu" in config["methodologie"] else ## "" if "otu" is not selected
        ["{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/MeBaPiNa_q2filter_uchime.report",
        "{tmp}02_analysis_results/03_otu_picking/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_counttax_q2kmermap.report"])

        ## ALIGN ##

        ([""] if not "align" in config["methodologie"] else ## "" if "align" is not selected
        ["{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/MeBaPiNa_filter_aligned.report",
        "{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_counttax_aligned.report"])

        ## K-MER ##

        ([""] if not "kmer" in config["methodologie"] else ## "" if "kmer" is not selected
        ["{tmp}01_processed_data/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_kmermap_filtered.report",
        "{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_retax_kmermap.report",
        "{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_counttax_kmermap.report"])

    )),
    tmp = config["experiments"]["tmp"],
    run = RUNS,
    barc = all_barcs,
    reference = config['reference']['source'],
    reftype = config['reference']['rank'] )
    ## return
    return input_list

rule all:
    input:
        input_barc
    params:
        config["experiments"]["tmp"] + "METADATA/ANALYSIS_PROGRESS_MANAGEMENT.csv"
    shell:
        "report_file={params}; "
        "if [[ ! -f ${{report_file}} ]]; then "
        "echo \"Sample name;File/directory;Completion date;Checksum;Performed by;Description\" > ${{report_file}}; fi; "
        "indiv_reports=( $(echo \"{input}\") ); "
        "for rprt in ${{indiv_reports[@]}}; do cat ${{rprt}} >> ${{report_file}}; done"




# def basecalls_per_barcode(wildcards):
#     ## get "pass" directory and trigger checkpoint (this way we can specify output inside the checkpoints output directory "pass" without direct rule association)
#     pass_dir = checkpoints.basecalling_raw.get(run=wildcards.run,tmp=wildcards.tmp).output[0]
#     ## get barcode for sample
#     sample_barcode = wildcards.barc
#     ## create file name for barcode
#     barc_input = pass_dir + "/" + sample_barcode
#     return barc_input
#
# rule all:
#     input:
#         expand(list(filter(None,[
#
#         ## BASECALL ##
#
#         ## general QC: all reads, including calibtation strads, intentional downsampling
#         "{tmp}02_analysis_results/01_basecalling/{run}/nanoplot/NanoStats.txt",
#         ## general QC: all reads, forced downsampling
#         "{tmp}02_analysis_results/01_basecalling/{run}/pycoqc/pycoQC_report.json",
#         ## per base QC: all reads, forced downsampling
#         "{tmp}02_analysis_results/01_basecalling/{run}/nanoqc/nanoQC.html",
#         ## read QC: all passed reads
#         "{tmp}02_analysis_results/01_basecalling/{run}/fastqc/stdin_fastqc.html",
#
#         ## barcode QC: per barcode
#         ("" if not BAC_KIT else ## "" if bac_kit is ""
#         "{tmp}02_analysis_results/01_basecalling/{run}/nanocomp/NanoStats.txt"),
#
#         ## TRIM AND FILTER ##
#
#         ## general QC: trimed and filtered barcoded reads, intentional downsampling
#         "{tmp}02_analysis_results/02_trimming_filtering/{run}/nanoplot/NanoStats.txt",
#         ## per base QC: trimed and filtered barcoded reads, forced downsampling
#         "{tmp}02_analysis_results/02_trimming_filtering/{run}/nanoqc/nanoQC.html",
#         ## read QC: trimed and filtered barcoded reads
#         "{tmp}02_analysis_results/02_trimming_filtering/{run}/fastqc/stdin_fastqc.html",
#
#         ## barcode QC: trimed and filtered barcoded reads
#         ("" if not BAC_KIT else ## "" if bac_kit is ""
#         "{tmp}02_analysis_results/02_trimming_filtering/{run}/nanocomp/NanoStats.txt"),
#
#         ## ALIGNMENT ##
#
#         ## general QC: per barcode, intentional downsampling
#         ("" if not "align" in config["methodologie"] else ## "" if "align" is not selected
#         "{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}_{reftype}/pycoqc.html"),
#         ## taxonomic composition
#         ("" if not "align" in config["methodologie"] else ## "" if "align" is not selected
#         "{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}_{reftype}/aligned.counttaxlist"),
#         ("" if not "align" in config["methodologie"] else ## "" if "align" is not selected
#         "{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}_{reftype}/krona.html"),
#
#         # K-MER MAPPING ##
#
#         ## taxonomic composition
#         ("" if not "kmer" in config["methodologie"] else ## "" if "kmer" is not selected
#         "{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/kmer.counttaxlist"),
#         ("" if not "kmer" in config["methodologie"] else ## "" if "kmer" is not selected
#         "{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/krona.html"),
#         ("" if not "kmer" in config["methodologie"] else ## "" if "kmer" is not selected
#         "{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/krona_bracken.html"),
#
#         ## OTU ##
#
#         ## clustered reads
#         ("" if not "otu" in config["methodologie"] else ## "" if "otu" is not selected
#         "{tmp}02_analysis_results/03_otu_picking/{run}/{barc}/{reference}/q2otupick/index.html"),
#         ## filtered reads
#         ("" if not "otu" in config["methodologie"] else ## "" if "otu" is not selected
#         "{tmp}02_analysis_results/03_otu_picking/{run}/{barc}/{reference}/q2filter/index.html"),
#         ## classified taxa
#         ("" if not "otu" in config["methodologie"] else ## "" if "otu" is not selected
#         "{tmp}02_analysis_results/03_otu_picking/{run}/{barc}/{reference}_{reftype}/kmer.counttaxlist"),
#         ("" if not "otu" in config["methodologie"] else ## "" if "otu" is not selected
#         "{tmp}02_analysis_results/03_otu_picking/{run}/{barc}/{reference}_{reftype}/krona.html"),
#
#         ## CALIBRATION STRAIN ##
#
#         ## calibration QC: only calinration strands
#         ("" if not LAM_DCS else ## "" if lam_DCS is False
#         "{tmp}02_analysis_results/01_basecalling/{run}_calibration_strands/nanoplot/NanoStats.txt"),
#
#         ## calibration QC: only calinration strands
#         ("" if not LAM_DCS else ## "" if lam_DCS is False
#         "{tmp}02_analysis_results/03_alignment/{run}_calibration_strands/lambda_nanoplot/NanoStats.txt"),
#         ("" if not LAM_DCS else ## "" if lam_DCS is False
#         "{tmp}02_analysis_results/03_alignment/{run}_calibration_strands/lambda_pycoqc/pycoQC_report.json"),
#
#         ## REPORT ##
#         "{tmp}00_raw_data/{run}_ANALYSIS_PROGRESS_MANAGEMENT.csv"
#
#         ])), tmp = config["experiments"]["tmp"], run = RUNS, barc = SAMPLES.keys(), reference = config['reference']['source'], reftype = config['reference']['rank'])
#
# ## REPORT ##
# ############
#
# rule report:
#     input:
#         "{tmp}00_raw_data/{run}/MeBaPiNa_moving_raw.report"
#     output:
#         temp("{tmp}00_raw_data/{run}_ANALYSIS_PROGRESS_MANAGEMENT.csv")
#     shell:
#         "{output}"
