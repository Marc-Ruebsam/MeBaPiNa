#############
## REPORTS ##
#############

## collection of all reports
def input_report(wildcards):
    from os import listdir
    ## get "pass" directory
    basecall_dir = checkpoints.basecall_raw.get(tmp=wildcards.tmp,run=wildcards.run).output[0]
    ## get barcode directory names within "pass" directory
    all_barcs = listdir(basecall_dir)
    ## retain only folders containing one of the selected barcodes (not unassigned)
    all_barcs = [barc for barc in all_barcs if barc in SAMPLES.keys()]
    ## create file names with barcodes
    input_list = [
        ## BASECALL ##
        ["{tmp}00_raw_data/{run}/MeBaPiNa_move_raw.report",
        "{tmp}00_raw_data/{run}/MeBaPiNa_basecall_raw.report"] +
        ## TRIM AND FILTER ##
        ["{tmp}01_processed_data/02_trimming_filtering/{run}/" + barc + "/MeBaPiNa_trim_basecalled.report",
        "{tmp}01_processed_data/02_trimming_filtering/{run}/" + barc + "/MeBaPiNa_filter_trimmed.report"] +
        ## OTU ##
        [x for x in
        ["{tmp}01_processed_data/03_otu_picking/{run}/" + barc + "/{reference}/MeBaPiNa_q2filter_uchime.report",
        "{tmp}02_analysis_results/03_otu_picking/{run}/" + barc + "/{reference}_{reftype}/MeBaPiNa_counttax_q2kmermap.report"]
        if "otu" in config["methodologie"]] +
        ## ALIGN ##
        [x for x in
        ["{tmp}01_processed_data/03_alignment/{run}/" + barc + "/{reference}/MeBaPiNa_filter_aligned.report",
        "{tmp}02_analysis_results/03_alignment/{run}/" + barc + "/{reference}_{reftype}/MeBaPiNa_counttax_aligned.report"]
        if "align" in config["methodologie"]] +
        ## K-MER ##
        [x for x in
        ["{tmp}01_processed_data/03_kmer_mapping/{run}/" + barc + "/{reference}_{reftype}/MeBaPiNa_kmermap_filtered.report",
        "{tmp}02_analysis_results/03_kmer_mapping/{run}/" + barc + "/{reference}_{reftype}/MeBaPiNa_retax_kmermap.report",
        "{tmp}02_analysis_results/03_kmer_mapping/{run}/" + barc + "/{reference}_{reftype}/MeBaPiNa_counttax_kmermap.report"]
        if "kmer" in config["methodologie"]]
        for barc in all_barcs]
    ## flatten list of lists
    input_list = [item for sublist in input_list for item in sublist]
    ## return
    return input_list

## rule for reports
rule all_report:
    input:
        input_report
    output:
        temp("{tmp}METADATA/{run}-{reference}-{reftype}.csv")
    shell:
        "report_file=$(echo \"{output}\" | sed 's#{wildcards.run}-{wildcards.reference}-{wildcards.reftype}#ANALYSIS_PROGRESS_MANAGEMENT#'); "
        "if [[ ! -f ${{report_file}} ]]; then "
        "echo \"Sample name;File/directory;Completion date;Checksum;Performed by;Description\" > ${{report_file}}; fi; "
        "indiv_reports=( $(echo \"{input}\") ); "
        "for rprt in ${{indiv_reports[@]}}; do cat ${{rprt}} >> ${{report_file}}; done; "
        "awk 'NR == 1; NR > 1 {{print $0 | \"sort -n | uniq\"}}' ${{report_file}} > {output}; " ## store unique lines in temporary output
        "cp {output} > ${{report_file}}" ## cp unique lines into output

## BASECALLING ##
#################

rule report_move_raw:
    input:
        "{tmp}00_raw_data/{run}/fast5"
    output:
        temp("{tmp}00_raw_data/{run}/MeBaPiNa_move_raw.report")
    params:
        " ".join(SAMPLES.values())
    shell:
        "all_IDs=( $(echo \"{params}\") ); "
        "for id in ${{all_IDs[@]}}; do echo "
        ## "Sample name;File/directory;Completion date;Checksum;Performed by;Description"
        "\"${{id}};{input};$(stat -c %y {input});NA;MeBaPiNa;General: moved raw fast5 files to path ready for analysis.\" "
        ">> {output}; done"

rule report_basecall_raw:
    input:
        "{tmp}01_processed_data/01_basecalling/{run}/sequencing_summary.txt",
        "{tmp}01_processed_data/01_basecalling/{run}/pass"
    output:
        temp("{tmp}00_raw_data/{run}/MeBaPiNa_basecall_raw.report")
    params:
        "all_IDs=( $(echo \"" + " ".join(SAMPLES.values()) + "\") ); "
    shell:
        "{params}"
        "for id in ${{all_IDs[@]}}; do echo "
        ## "Sample name;File/directory;Completion date;Checksum;Performed by;Description"
        "\"${{id}};{input};$(stat -c %y {input[0]});NA;MeBaPiNa;General: basecalled and demultiplexed raw fast5 files into fastq.\" "
        ">> {output}; done"

## TRIM AND FILTER ##
#####################

rule report_trim_basecalled:
    input:
        "{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/trimmed.fastq"
    output:
        temp("{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/MeBaPiNa_trim_basecalled.report")
    params:
        "all_IDs=( $(echo \"" + " ".join(SAMPLES.values()) + "\") ); ",
        "all_barcs=( $(echo \"" + " ".join(SAMPLES.keys()) + "\") ); "
    shell:
        "{params}"
        "id=${{all_IDs[$(for i in \"${{!all_barcs[@]}}\"; do "
        "if [[ \"${{all_barcs[$i]}}\" = \"{wildcards.barc}\" ]]; then echo $i; fi; "
        "done)]}}; echo "
        ## "Sample name;File/directory;Completion date;Checksum;Performed by;Description"
        "\"${{id}};{input};$(stat -c %y {input});NA;MeBaPiNa;General: trimmed adapters and barcodes from reads and demultiplexed a second time.\" "
        ">> {output}"

rule report_filter_trimmed:
    input:
        "{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/filtered.fastq"
    output:
        temp("{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/MeBaPiNa_filter_trimmed.report")
    params:
        "all_IDs=( $(echo \"" + " ".join(SAMPLES.values()) + "\") ); ",
        "all_barcs=( $(echo \"" + " ".join(SAMPLES.keys()) + "\") ); "
    shell:
        "{params}"
        "id=${{all_IDs[$(for i in \"${{!all_barcs[@]}}\"; do "
        "if [[ \"${{all_barcs[$i]}}\" = \"{wildcards.barc}\" ]]; then echo $i; fi; "
        "done)]}}; echo "
        ## "Sample name;File/directory;Completion date;Checksum;Performed by;Description"
        "\"${{id}};{input};$(stat -c %y {input});NA;MeBaPiNa;General: length and quality filtered reads.\" "
        ">> {output}"

## OTU ##
#########

rule report_q2filter_uchime:
    input:
        nochimtable="{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/filt_ftable.qza",
        nochimseq="{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/filt_centseq.qza"
    output:
        temp("{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/MeBaPiNa_q2filter_uchime.report")
    params:
        "all_IDs=( $(echo \"" + " ".join(SAMPLES.values()) + "\") ); ",
        "all_barcs=( $(echo \"" + " ".join(SAMPLES.keys()) + "\") ); "
    shell:
        "{params}"
        "id=${{all_IDs[$(for i in \"${{!all_barcs[@]}}\"; do "
        "if [[ \"${{all_barcs[$i]}}\" = \"{wildcards.barc}\" ]]; then echo $i; fi; "
        "done)]}}; echo "
        ## "Sample name;File/directory;Completion date;Checksum;Performed by;Description"
        "\"${{id}};{input};$(stat -c %y {input[0]});NA;MeBaPiNa;OTU: dereplication, open-reference clustering, chimera removal and filtering.\" "
        ">> {output}"

rule report_counttax_q2kmermap:
    input:
        counttaxlist="{tmp}02_analysis_results/03_otu_picking/{run}/{barc}/{reference}_{reftype}/kmer.counttaxlist"
    output:
        temp("{tmp}02_analysis_results/03_otu_picking/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_counttax_q2kmermap.report")
    params:
        "all_IDs=( $(echo \"" + " ".join(SAMPLES.values()) + "\") ); ",
        "all_barcs=( $(echo \"" + " ".join(SAMPLES.keys()) + "\") ); "
    shell:
        "{params}"
        "id=${{all_IDs[$(for i in \"${{!all_barcs[@]}}\"; do "
        "if [[ \"${{all_barcs[$i]}}\" = \"{wildcards.barc}\" ]]; then echo $i; fi; "
        "done)]}}; echo "
        ## "Sample name;File/directory;Completion date;Checksum;Performed by;Description"
        "\"${{id}};{input};$(stat -c %y {input[0]});NA;MeBaPiNa;OTU: taxonomic classification and file conversion.\" "
        ">> {output}"

## ALIGN ##
###########

rule report_filter_aligned:
    input:
        bam="{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/filteredsorted.bam",
        bai="{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/filteredsorted.bam.bai"
    output:
        temp("{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/MeBaPiNa_filter_aligned.report")
    params:
        "all_IDs=( $(echo \"" + " ".join(SAMPLES.values()) + "\") ); ",
        "all_barcs=( $(echo \"" + " ".join(SAMPLES.keys()) + "\") ); "
    shell:
        "{params}"
        "id=${{all_IDs[$(for i in \"${{!all_barcs[@]}}\"; do "
        "if [[ \"${{all_barcs[$i]}}\" = \"{wildcards.barc}\" ]]; then echo $i; fi; "
        "done)]}}; echo "
        ## "Sample name;File/directory;Completion date;Checksum;Performed by;Description"
        "\"${{id}};{input};$(stat -c %y {input[0]});NA;MeBaPiNa;Alignment: alignment, extraction of uniquely aligned reads (not fully filtered).\" "
        ">> {output}"

rule report_counttax_aligned:
    input:
        counttaxlist="{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}_{reftype}/aligned.counttaxlist"
    output:
        temp("{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_counttax_aligned.report")
    params:
        "all_IDs=( $(echo \"" + " ".join(SAMPLES.values()) + "\") ); ",
        "all_barcs=( $(echo \"" + " ".join(SAMPLES.keys()) + "\") ); "
    shell:
        "{params}"
        "id=${{all_IDs[$(for i in \"${{!all_barcs[@]}}\"; do "
        "if [[ \"${{all_barcs[$i]}}\" = \"{wildcards.barc}\" ]]; then echo $i; fi; "
        "done)]}}; echo "
        ## "Sample name;File/directory;Completion date;Checksum;Performed by;Description"
        "\"${{id}};{input};$(stat -c %y {input[0]});NA;MeBaPiNa;Alignment: filtering and taxonomic classification.\" "
        ">> {output}"

## K-MER MAPPING ##
###################

rule report_kmermap_filtered:
    input:
        report="{tmp}01_processed_data/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/filtered.kreport2",
        output="{tmp}01_processed_data/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/filtered.kraken2"
    output:
        temp("{tmp}01_processed_data/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_kmermap_filtered.report")
    params:
        "all_IDs=( $(echo \"" + " ".join(SAMPLES.values()) + "\") ); ",
        "all_barcs=( $(echo \"" + " ".join(SAMPLES.keys()) + "\") ); "
    shell:
        "{params}"
        "id=${{all_IDs[$(for i in \"${{!all_barcs[@]}}\"; do "
        "if [[ \"${{all_barcs[$i]}}\" = \"{wildcards.barc}\" ]]; then echo $i; fi; "
        "done)]}}; echo "
        ## "Sample name;File/directory;Completion date;Checksum;Performed by;Description"
        "\"${{id}};{input};$(stat -c %y {input[0]});NA;MeBaPiNa;K-mer: taxonomic classification.\" "
        ">> {output}"

rule report_retax_kmermap:
    input:
        "{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/{reftype}.bracken",
        "{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/{reftype}.kreport2",
    output:
        temp("{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_retax_kmermap.report")
    params:
        "all_IDs=( $(echo \"" + " ".join(SAMPLES.values()) + "\") ); ",
        "all_barcs=( $(echo \"" + " ".join(SAMPLES.keys()) + "\") ); "
    shell:
        "{params}"
        "id=${{all_IDs[$(for i in \"${{!all_barcs[@]}}\"; do "
        "if [[ \"${{all_barcs[$i]}}\" = \"{wildcards.barc}\" ]]; then echo $i; fi; "
        "done)]}}; echo "
        ## "Sample name;File/directory;Completion date;Checksum;Performed by;Description"
        "\"${{id}};{input};$(stat -c %y {input[0]});NA;MeBaPiNa;K-mer: abundance reestimation.\" "
        ">> {output}"

rule report_counttax_kmermap:
    input:
        counttaxlist="{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/kmer.counttaxlist"
    output:
        temp("{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_counttax_kmermap.report")
    params:
        "all_IDs=( $(echo \"" + " ".join(SAMPLES.values()) + "\") ); ",
        "all_barcs=( $(echo \"" + " ".join(SAMPLES.keys()) + "\") ); "
    shell:
        "{params}"
        "id=${{all_IDs[$(for i in \"${{!all_barcs[@]}}\"; do "
        "if [[ \"${{all_barcs[$i]}}\" = \"{wildcards.barc}\" ]]; then echo $i; fi; "
        "done)]}}; echo "
        ## "Sample name;File/directory;Completion date;Checksum;Performed by;Description"
        "\"${{id}};{input};$(stat -c %y {input[0]});NA;MeBaPiNa;K-mer: fle conversion.\" "
        ">> {output}"

## CALIBRATION STRAIN ##
########################
