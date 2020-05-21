#############
## REPORTS ##
#############

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
        "\"${{id}};{input};$(stat -c %y {input});NA;MeBaPiNa;Moved raw fast5 files to path ready for analysis.\" "
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
        "\"${{id}};{input};$(stat -c %y {input[0]});NA;MeBaPiNa;Basecalled and demultiplexed raw fast5 files into fastq.\" "
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
        "\"${{id}};{input};$(stat -c %y {input});NA;MeBaPiNa;Trimmed adapters and barcodes from reads and demultiplexed a second time.\" "
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
        "\"${{id}};{input};$(stat -c %y {input});NA;MeBaPiNa;Length and quality filtered reads.\" "
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
