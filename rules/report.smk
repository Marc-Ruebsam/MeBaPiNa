##########
## LOGS ##
##########

## basecalling

rule moving_raw_basecalling:
    input:
        "{tmp}00_raw_data/{run}/fast5"
    output:
        "test.test"
    params:
        SAMPELS
    shell:
        "{params}"



















############################
## CHECKPOINT AGGREGATION ##
############################

# def input_aggregate(wildcards):
#     from os import listdir
#     ## "pass" directory
#     basecall_dir = checkpoints.guppy.get(run=wildcards.run).output[0]
#     ## directory names within "pass" directory
#     barcs = listdir(basecall_dir)
#     ## retain only strings containing "barcode"
#     barcs = [i for i in barcs if "barcode" in i]
#     ## create file names with barcodes
#     barc_input = expand("{tmp}02_analysis_results/03_alignment/{run}/{barc}_pycoqc/pycoQC_report.json",
#         run=wildcards.run,
#         barc=barcs)
#     barc_input.sort()
#     return barc_input
#
# rule aggr_align_barc:
#     input:
#         input_aggregate ## if it keeps returning only the pass folder, some but not all of the output and input of the guppy rule exists. Create it, even if its only an empty folder to fix the problem.
#     output:
#         "{tmp}02_analysis_results/03_alignment/{run}/MeBaPiNa_barcode_aggregation.txt"
#     shell:
#         "echo {input} > {output}"
