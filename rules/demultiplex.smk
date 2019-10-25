# rule porechop:
#     input:
#         "01_processeddata/{run}/pass"
#     output:
#         directory("01_processeddata/{run}/basecall_demultiplex")
#     log:
#         "01_processeddata/{run}/pass_demultiplexed/MeBaPiNa_porechop.log" ##  > {log} 2>&1 (TLDR both to log) redirects stdout to file and stderr to stdout which is redirected to a file
#     benchmark:
#         "01_processeddata/{run}/pass_demultiplexed/MeBaPiNa_porechop.benchmark.tsv"
#     conda:
#         "../envs/porechop.yml"
#     threads:
#         config["machine"]["cpu"]
#     params:
#         "--verbosity 2" ## or nothing to log
#     shell:
#         "porechop --threads {threads} {params} --input {input} --barcode_dir {output} > {log} 2>&1" 
# 
# rule guppy_demulti:
#     input:
#         "01_processeddata/{run}/basecall/pass"
#     output:
#         "01_processeddata/{run}/basecall_demultiplex/barcoding_summary.txt"
#     log:
#         "01_processeddata/{run}/basecall_demultiplex/MeBaPiNa_guppy_demulti.log" ##  > {log} 2>&1 (TLDR both to log) redirects stdout to file and stderr to stdout which is redirected to a file
#     benchmark:
#         "01_processeddata/{run}/basecall_demultiplex/MeBaPiNa_guppy_demulti.benchmark.tsv"
#     version:
#         subprocess.check_output("guppy_barcoder --version | awk '{print $NF}'", shell=True)
#     threads:
#         config["machine"]["cpu"]
#     params:
#         ("--barcode_kits " + config["guppy"]["bac_kit"]), ## includes demultiplexing
#         ("--device cuda:all:100%" if config["machine"]["gpu"] else "")
#         # --require_barcodes_both_ends
#         # --detect_mid_strand_barcodes ## slower but mid barcodes
#         # --trim_barcodes ## trim barcodes
#         # --num_extra_bases_trim arg ## trim extra bases
#     shell:
#         "guppy_barcoder --worker_threads {threads} {params} "
#         "--save_path {output} "
#         "--input_path {input} > {log} 2>&1"
# 
# ruleorder: guppy_demulti > porechop
