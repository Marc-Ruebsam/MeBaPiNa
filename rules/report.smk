#############
## REPORTS ##
#############

## COPY OUTPUT ##
#################

## BASECALL ##

rule copy_basecall_output:
    input:
        "{tmp}02_analysis_results/01_basecalling/{run}/nanoplot/NanoPlot-report.html", ## general QC: all reads, including calibtation strads, intentional downsampling
        "{tmp}02_analysis_results/01_basecalling/{run}/nanoplot/NanoStats.txt", ## general QC: all reads, including calibtation strads, intentional downsampling
        "{tmp}02_analysis_results/01_basecalling/{run}/pycoqc/pycoQC_report.html", ## general QC: all reads, forced downsampling
        "{tmp}02_analysis_results/01_basecalling/{run}/pycoqc/pycoQC_report.json", ## general QC: all reads, forced downsampling
        "{tmp}02_analysis_results/01_basecalling/{run}/nanocomp/NanoComp-report.html", ## barcode QC: per barcode
        "{tmp}02_analysis_results/01_basecalling/{run}/nanocomp/NanoStats.txt", ## barcode QC: per barcode
        "{tmp}02_analysis_results/01_basecalling/{run}/nanoqc/nanoQC.html", ## per base QC: all reads, forced downsampling
        "{tmp}02_analysis_results/01_basecalling/{run}/fastqc/stdin_fastqc.html", ## read QC: all passed reads
        "{tmp}00_raw_data/{run}/MeBaPiNa_move_raw.report", ## REPORT
        "{tmp}00_raw_data/{run}/MeBaPiNa_basecall_raw_seqsum.report", ## REPORT
        "{tmp}00_raw_data/{run}/MeBaPiNa_basecall_raw_pass.report" ## REPORT
    output:
        "{tmp}03_report/{timepoint}/{sample}/{run}-{barc}/01_basecalling-nanoplot-NanoPlot-report.html",
        "{tmp}03_report/{timepoint}/{sample}/{run}-{barc}/01_basecalling-nanoplot-NanoStats.txt",
        "{tmp}03_report/{timepoint}/{sample}/{run}-{barc}/01_basecalling-pycoqc-pycoQC_report.html",
        "{tmp}03_report/{timepoint}/{sample}/{run}-{barc}/01_basecalling-pycoqc-pycoQC_report.json",
        "{tmp}03_report/{timepoint}/{sample}/{run}-{barc}/01_basecalling-nanocomp-NanoComp-report.html",
        "{tmp}03_report/{timepoint}/{sample}/{run}-{barc}/01_basecalling-nanocomp-NanoStats.txt",
        "{tmp}03_report/{timepoint}/{sample}/{run}-{barc}/01_basecalling-nanoqc-nanoQC.html",
        "{tmp}03_report/{timepoint}/{sample}/{run}-{barc}/01_basecalling-fastqc-stdin_fastqc.html"
    shell:
        "all_input=( $(echo \"{input}\") ); "
        "all_output=( $(echo \"{output}\") ); "
        "for idx in \"${{!all_output[@]}}\"; do " ## loop over index of output array (because of dummy report files)
        "cp \"${{all_input[$idx]}}\" \"${{all_output[$idx]}}\"; done" ## copy input at index idx to output at index idx

## TRIM AND FILTER ##

rule copy_trim_filter_output:
    input:
        "{tmp}02_analysis_results/02_trimming_filtering/{run}/nanoplot/NanoPlot-report.html", ## general QC: trimed and filtered barcoded reads, intentional downsampling
        "{tmp}02_analysis_results/02_trimming_filtering/{run}/nanoplot/NanoStats.txt", ## general QC: trimed and filtered barcoded reads, intentional downsampling
        "{tmp}02_analysis_results/02_trimming_filtering/{run}/nanocomp/NanoComp-report.html", ## barcode QC: trimed and filtered barcoded reads
        "{tmp}02_analysis_results/02_trimming_filtering/{run}/nanocomp/NanoStats.txt", ## barcode QC: trimed and filtered barcoded reads
        "{tmp}02_analysis_results/02_trimming_filtering/{run}/nanoqc/nanoQC.html", ## per base QC: trimed and filtered barcoded reads, forced downsampling
        "{tmp}02_analysis_results/02_trimming_filtering/{run}/fastqc/stdin_fastqc.html", ## read QC: trimed and filtered barcoded reads
        "{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/MeBaPiNa_trim_basecalled.report", ## REPORT
        "{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/MeBaPiNa_filter_trimmed.report" ## REPORT
    output:
        "{tmp}03_report/{timepoint}/{sample}/{run}-{barc}/02_trimming_filtering-nanoplot-NanoPlot-report.html",
        "{tmp}03_report/{timepoint}/{sample}/{run}-{barc}/02_trimming_filtering-nanoplot-NanoStats.txt",
        "{tmp}03_report/{timepoint}/{sample}/{run}-{barc}/02_trimming_filtering-nanocomp-NanoComp-report.html",
        "{tmp}03_report/{timepoint}/{sample}/{run}-{barc}/02_trimming_filtering-nanocomp-NanoStats.txt",
        "{tmp}03_report/{timepoint}/{sample}/{run}-{barc}/02_trimming_filtering-nanoqc-nanoQC.html",
        "{tmp}03_report/{timepoint}/{sample}/{run}-{barc}/02_trimming_filtering-fastqc-stdin_fastqc.html"
    shell:
        "all_input=( $(echo \"{input}\") ); "
        "all_output=( $(echo \"{output}\") ); "
        "for idx in \"${{!all_output[@]}}\"; do "
        "cp \"${{all_input[$idx]}}\" \"${{all_output[$idx]}}\"; done"

## OTU ##

rule copy_otu_output_I: ## because of different wildcards in output
    input:
        ## input to copy
        "{tmp}02_analysis_results/03_otu_picking/{run}/{barc}/{reference}/q2otupick/index.html", ## clustered reads
        "{tmp}02_analysis_results/03_otu_picking/{run}/{barc}/{reference}/q2filter/index.html", ## filtered reads
        ## dummy depencencies
        "{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/MeBaPiNa_q2filter_uchime_ftable.report", ## REPORT
        "{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/MeBaPiNa_q2filter_uchime_centseq.report" ## REPORT
    output:
        "{tmp}03_report/{timepoint}/{sample}/{run}-{barc}/03_otu_picking-{reference}-q2otupick-index.html",
        "{tmp}03_report/{timepoint}/{sample}/{run}-{barc}/03_otu_picking-{reference}-q2filter-index.html"
    shell:
        "all_input=( $(echo \"{input}\") ); "
        "all_output=( $(echo \"{output}\") ); "
        "for idx in \"${{!all_output[@]}}\"; do "
        "cp \"${{all_input[$idx]}}\" \"${{all_output[$idx]}}\"; done"

rule copy_otu_output_II: ## because of different wildcards in output
    input:
        ## input to copy
        "{tmp}02_analysis_results/03_otu_picking/{run}/{barc}/{reference}_{reftype}/krona.html", ## classified taxa
        "{tmp}02_analysis_results/03_otu_picking/{run}/{barc}/{reference}_{reftype}/kmer.counttaxlist", ## taxonomic classifications
        ## dummy depencencies
        "{tmp}02_analysis_results/03_otu_picking/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_kmermap_q2rereplicate.report", ## REPORT
        "{tmp}02_analysis_results/03_otu_picking/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_counttax_q2kmermap.report" ## REPORT
    output:
        "{tmp}03_report/{timepoint}/{sample}/{run}-{barc}/03_otu_picking-{reference}_{reftype}-krona.html",
        "{tmp}03_report/{timepoint}/{sample}/{run}-{barc}/03_otu_picking-{reference}_{reftype}-kmer.counttaxlist"
    shell:
        "all_input=( $(echo \"{input}\") ); "
        "all_output=( $(echo \"{output}\") ); "
        "for idx in \"${{!all_output[@]}}\"; do "
        "cp \"${{all_input[$idx]}}\" \"${{all_output[$idx]}}\"; done"

## ALIGNMENT ##

rule copy_align_output:
    input:
        ## input to copy
        "{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}_{reftype}/pycoqc.html", ## per barcode, intentional downsampling
        "{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}_{reftype}/pycoqc.json", ## per barcode, intentional downsampling
        "{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}_{reftype}/krona.html", ## taxonomic classification
        "{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}_{reftype}/aligned.counttaxlist", ## taxonomic classification
        ## dummy depencencies
        "{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/MeBaPiNa_filter_aligned.report", ## REPORT
        "{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_counttax_aligned.report" ## REPORT
    output:
        "{tmp}03_report/{timepoint}/{sample}/{run}-{barc}/03_alignment/{reference}_{reftype}-pycoqc.html",
        "{tmp}03_report/{timepoint}/{sample}/{run}-{barc}/03_alignment/{reference}_{reftype}-pycoqc.json",
        "{tmp}03_report/{timepoint}/{sample}/{run}-{barc}/03_alignment/{reference}_{reftype}-krona.html",
        "{tmp}03_report/{timepoint}/{sample}/{run}-{barc}/03_alignment/{reference}_{reftype}-aligned.counttaxlist"
    shell:
        "all_input=( $(echo \"{input}\") ); "
        "all_output=( $(echo \"{output}\") ); "
        "for idx in \"${{!all_output[@]}}\"; do "
        "cp \"${{all_input[$idx]}}\" \"${{all_output[$idx]}}\"; done"

## K-MER MAPPING ##

rule copy_kmer_output:
    input:
        ## input to copy
        "{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/krona.html", ## taxonomic composition
        "{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/krona_bracken.html", ## taxonomic composition after reestimation
        "{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/kmer.counttaxlist", ## taxonomic classification
        ## dummy depencencies
        "{tmp}01_processed_data/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_kmermap_filtered.report", ## REPORT
        "{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_retax_kmermap.report", ## REPORT
        "{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_counttax_kmermap.report" ## REPORT
    output:
        "{tmp}03_report/{timepoint}/{sample}/{run}-{barc}/03_kmer_mapping/{reference}_{reftype}-krona.html",
        "{tmp}03_report/{timepoint}/{sample}/{run}-{barc}/03_kmer_mapping/{reference}_{reftype}-krona_bracken.html",
        "{tmp}03_report/{timepoint}/{sample}/{run}-{barc}/03_kmer_mapping/{reference}_{reftype}-kmer.counttaxlist"
    shell:
        "all_input=( $(echo \"{input}\") ); "
        "all_output=( $(echo \"{output}\") ); "
        "for idx in \"${{!all_output[@]}}\"; do "
        "cp \"${{all_input[$idx]}}\" \"${{all_output[$idx]}}\"; done"


## CREATE REPORT FILE ##
########################

rule report_initiate:
    output:
        LOGS
    shell:
        "if [[ ! -f \"{output}\" ]]; then "
        "echo \"Sample name;File/Directory;Action;Date;Checksum;Who;Description\""
        " > \"{output}\"; else touch \"{output}\"; fi"

## OUTPUT CREATION REPORTS ##
#############################

## BASECALLING ##

rule report_move_raw:
    input:
        thing="{tmp}00_raw_data/{run}/fast5"
    output:
        dummy=temp("{tmp}00_raw_data/{run}/MeBaPiNa_move_raw.report") ## for dependencies
    params:
        "report_log=" + LOGS + "; ", ## text directly reported here
        "all_IDs=( $(echo \"" + " ".join(SAMPLES.values()) + "\") ); " ## will be reported for all samples in this run
    shell:
        "{params}"
        "md5checksum=$(find \"{input.thing}\" -type f -exec md5sum {{}} \\; | sort -k 2 | awk '{{print $1}}' | md5sum | awk '{{print $1}}'); " ## md5 of md5s
        "for id in ${{all_IDs[@]}}; do echo "
        "\"${{id}};{input.thing};Moved input;$(stat -c '%.19z' {input.thing});${{md5checksum}};MeBaPiNa;"
        "General: raw fast5 files moved to directory ready for analysis. Checksum was calculated as: find -type f -exec md5sum {{}} \\; | sort -k 2 | awk '{{print \$1}}' | md5sum \" "
        ">> \"${{report_log}}\"; done; "
        "touch {output.dummy}" ## for dependencies

rule report_basecall_raw_seqsum:
    input:
        thing="{tmp}01_processed_data/01_basecalling/{run}/sequencing_summary.txt"
    output:
        dummy=temp("{tmp}00_raw_data/{run}/MeBaPiNa_basecall_raw_seqsum.report") ## for dependencies
    params:
        "report_log=" + LOGS + "; ", ## text directly reported here
        "all_IDs=( $(echo \"" + " ".join(SAMPLES.values()) + "\") ); " ## will be reported for all samples in this run
    shell:
        "{params}"
        "md5checksum=$(md5sum \"{input.thing}\" | awk '{{print $1}}'); "
        "for id in ${{all_IDs[@]}}; do echo "
        "\"${{id}};{input.thing};Created output;$(stat -c '%.19z' {input.thing});${{md5checksum}};MeBaPiNa;"
        "General: basecalling and demultiplexing log per read.\" "
        ">> \"${{report_log}}\"; done; "
        "touch {output.dummy}" ## for dependencies

rule report_basecall_raw_pass:
    input:
        thing="{tmp}01_processed_data/01_basecalling/{run}/pass"
    output:
        dummy=temp("{tmp}00_raw_data/{run}/MeBaPiNa_basecall_raw_pass.report")
    params:
        "report_log=" + LOGS + "; ",
        "all_IDs=( $(echo \"" + " ".join(SAMPLES.values()) + "\") ); "
    shell:
        "{params}"
        "md5checksum=$(find \"{input.thing}\" -type f -exec md5sum {{}} \\; | sort -k 2 | awk '{{print $1}}' | md5sum | awk '{{print $1}}'); "
        "for id in ${{all_IDs[@]}}; do echo "
        "\"${{id}};{input.thing};Created output;$(stat -c '%.19z' {input.thing});${{md5checksum}};MeBaPiNa;"
        "General: basecalled and demultiplexed fastq files. Checksum was calculated as: find -type f -exec md5sum {{}} \\; | sort -k 2 | awk '{{print \$1}}' | md5sum \" "
        ">> \"${{report_log}}\"; done; "
        "touch {output.dummy}"

## TRIM AND FILTER ##

rule report_trim_basecalled:
    input:
        thing="{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/trimmed.fastq"
    output:
        dummy=temp("{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/MeBaPiNa_trim_basecalled.report")
    params:
        "report_log=" + LOGS + "; ",
        "all_IDs=( $(echo \"" + " ".join(SAMPLES.values()) + "\") ); ",
        "all_barcs=( $(echo \"" + " ".join(SAMPLES.keys()) + "\") ); "
    shell:
        "{params}"
        "id=${{all_IDs[$(for idx in \"${{!all_barcs[@]}}\"; do " ## have to find sample name (id) to barcode and are not allowed wildcards as variable keys
        "if [[ \"${{all_barcs[$idx]}}\" = \"{wildcards.barc}\" ]]; then echo $idx; fi; done)]}}; "
        "md5checksum=$(md5sum \"{input.thing}\" | awk '{{print $1}}'); "
        "echo \"${{id}};{input.thing};Created output;$(stat -c '%.19z' {input.thing});${{md5checksum}};MeBaPiNa;"
        "General: trimmed adapters and barcodes from reads and demultiplexed a second time.\" "
        ">> \"${{report_log}}\"; "
        "touch {output.dummy}"

rule report_filter_trimmed:
    input:
        thing="{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/filtered.fastq"
    output:
        dummy=temp("{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/MeBaPiNa_filter_trimmed.report")
    params:
        "report_log=" + LOGS + "; ",
        "all_IDs=( $(echo \"" + " ".join(SAMPLES.values()) + "\") ); ",
        "all_barcs=( $(echo \"" + " ".join(SAMPLES.keys()) + "\") ); "
    shell:
        "{params}"
        "id=${{all_IDs[$(for idx in \"${{!all_barcs[@]}}\"; do "
        "if [[ \"${{all_barcs[$idx]}}\" = \"{wildcards.barc}\" ]]; then echo $idx; fi; done)]}}; "
        "md5checksum=$(md5sum \"{input.thing}\" | awk '{{print $1}}'); "
        "echo \"${{id}};{input.thing};Created output;$(stat -c '%.19z' {input.thing});${{md5checksum}};MeBaPiNa;"
        "General: length and quality filtered reads.\" "
        ">> \"${{report_log}}\"; "
        "touch {output.dummy}"

## OTU ##

rule report_q2filter_uchime_ftable:
    input:
        thing="{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/filt_ftable.qza"
    output:
        dummy=temp("{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/MeBaPiNa_q2filter_uchime_ftable.report")
    params:
        "report_log=" + LOGS + "; ",
        "all_IDs=( $(echo \"" + " ".join(SAMPLES.values()) + "\") ); ",
        "all_barcs=( $(echo \"" + " ".join(SAMPLES.keys()) + "\") ); "
    shell:
        "{params}"
        "id=${{all_IDs[$(for idx in \"${{!all_barcs[@]}}\"; do "
        "if [[ \"${{all_barcs[$idx]}}\" = \"{wildcards.barc}\" ]]; then echo $idx; fi; done)]}}; "
        "md5checksum=$(md5sum \"{input.thing}\" | awk '{{print $1}}'); "
        "echo \"${{id}};{input.thing};Created output;$(stat -c '%.19z' {input.thing});${{md5checksum}};MeBaPiNa;"
        "OTU: count table of features after: dereplication, open-reference clustering, chimera removal and filtering.\" "
        ">> \"${{report_log}}\"; "
        "touch {output.dummy}"

rule report_q2filter_uchime_centseq:
    input:
        thing="{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/filt_centseq.qza"
    output:
        dummy=temp("{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/MeBaPiNa_q2filter_uchime_centseq.report")
    params:
        "report_log=" + LOGS + "; ",
        "all_IDs=( $(echo \"" + " ".join(SAMPLES.values()) + "\") ); ",
        "all_barcs=( $(echo \"" + " ".join(SAMPLES.keys()) + "\") ); "
    shell:
        "{params}"
        "id=${{all_IDs[$(for idx in \"${{!all_barcs[@]}}\"; do "
        "if [[ \"${{all_barcs[$idx]}}\" = \"{wildcards.barc}\" ]]; then echo $idx; fi; done)]}}; "
        "md5checksum=$(md5sum \"{input.thing}\" | awk '{{print $1}}'); "
        "echo \"${{id}};{input.thing};Created output;$(stat -c '%.19z' {input.thing});${{md5checksum}};MeBaPiNa;"
        "OTU: feature clusters central sequences after: dereplication, open-reference clustering, chimera removal and filtering.\" "
        ">> \"${{report_log}}\"; "
        "touch {output.dummy}"

rule report_kmermap_q2rereplicate:
    input:
        thing="{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}_{reftype}/filtered.kreport2"
    output:
        dummy=temp("{tmp}02_analysis_results/03_otu_picking/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_kmermap_q2rereplicate.report")
    params:
        "report_log=" + LOGS + "; ",
        "all_IDs=( $(echo \"" + " ".join(SAMPLES.values()) + "\") ); ",
        "all_barcs=( $(echo \"" + " ".join(SAMPLES.keys()) + "\") ); "
    shell:
        "{params}"
        "id=${{all_IDs[$(for idx in \"${{!all_barcs[@]}}\"; do "
        "if [[ \"${{all_barcs[$idx]}}\" = \"{wildcards.barc}\" ]]; then echo $idx; fi; done)]}}; "
        "md5checksum=$(md5sum \"{input.thing}\" | awk '{{print $1}}'); "
        "echo \"${{id}};{input.thing};Created output;$(stat -c '%.19z' {input.thing});${{md5checksum}};MeBaPiNa;"
        "OTU: taxonomic classification.\" "
        ">> \"${{report_log}}\"; "
        "touch {output.dummy}"

rule report_counttax_q2kmermap:
    input:
        thing="{tmp}02_analysis_results/03_otu_picking/{run}/{barc}/{reference}_{reftype}/kmer.counttaxlist"
    output:
        dummy=temp("{tmp}02_analysis_results/03_otu_picking/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_counttax_q2kmermap.report")
    params:
        "report_log=" + LOGS + "; ",
        "all_IDs=( $(echo \"" + " ".join(SAMPLES.values()) + "\") ); ",
        "all_barcs=( $(echo \"" + " ".join(SAMPLES.keys()) + "\") ); "
    shell:
        "{params}"
        "id=${{all_IDs[$(for idx in \"${{!all_barcs[@]}}\"; do "
        "if [[ \"${{all_barcs[$idx]}}\" = \"{wildcards.barc}\" ]]; then echo $idx; fi; done)]}}; "
        "md5checksum=$(md5sum \"{input.thing}\" | awk '{{print $1}}'); "
        "echo \"${{id}};{input.thing};Created output;$(stat -c '%.19z' {input.thing});${{md5checksum}};MeBaPiNa;"
        "OTU: file conversion.\" "
        ">> \"${{report_log}}\"; "
        "touch {output.dummy}"

## ALIGN ##

rule report_filter_aligned:
    input:
        thing="{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/filteredsorted.bam"
    output:
        dummy=temp("{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/MeBaPiNa_filter_aligned.report")
    params:
        "report_log=" + LOGS + "; ",
        "all_IDs=( $(echo \"" + " ".join(SAMPLES.values()) + "\") ); ",
        "all_barcs=( $(echo \"" + " ".join(SAMPLES.keys()) + "\") ); "
    shell:
        "{params}"
        "id=${{all_IDs[$(for idx in \"${{!all_barcs[@]}}\"; do "
        "if [[ \"${{all_barcs[$idx]}}\" = \"{wildcards.barc}\" ]]; then echo $idx; fi; done)]}}; "
        "md5checksum=$(md5sum \"{input.thing}\" | awk '{{print $1}}'); "
        "echo \"${{id}};{input.thing};Created output;$(stat -c '%.19z' {input.thing});${{md5checksum}};MeBaPiNa;"
        "Alignment: alignment, extraction of uniquely aligned reads (not fully filtered).\" "
        ">> \"${{report_log}}\"; "
        "touch {output.dummy}"

rule report_counttax_aligned:
    input:
        thing="{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}_{reftype}/aligned.counttaxlist"
    output:
        dummy=temp("{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_counttax_aligned.report")
    params:
        "report_log=" + LOGS + "; ",
        "all_IDs=( $(echo \"" + " ".join(SAMPLES.values()) + "\") ); ",
        "all_barcs=( $(echo \"" + " ".join(SAMPLES.keys()) + "\") ); "
    shell:
        "{params}"
        "id=${{all_IDs[$(for idx in \"${{!all_barcs[@]}}\"; do "
        "if [[ \"${{all_barcs[$idx]}}\" = \"{wildcards.barc}\" ]]; then echo $idx; fi; done)]}}; "
        "md5checksum=$(md5sum \"{input.thing}\" | awk '{{print $1}}'); "
        "echo \"${{id}};{input.thing};Created output;$(stat -c '%.19z' {input.thing});${{md5checksum}};MeBaPiNa;"
        "Alignment: filtering, taxonomic classification and file conversion.\" "
        ">> \"${{report_log}}\"; "
        "touch {output.dummy}"

## K-MER MAPPING ##

rule report_kmermap_filtered:
    input:
        thing="{tmp}01_processed_data/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/filtered.kreport2"
    output:
        dummy=temp("{tmp}01_processed_data/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_kmermap_filtered.report")
    params:
        "report_log=" + LOGS + "; ",
        "all_IDs=( $(echo \"" + " ".join(SAMPLES.values()) + "\") ); ",
        "all_barcs=( $(echo \"" + " ".join(SAMPLES.keys()) + "\") ); "
    shell:
        "{params}"
        "id=${{all_IDs[$(for idx in \"${{!all_barcs[@]}}\"; do "
        "if [[ \"${{all_barcs[$idx]}}\" = \"{wildcards.barc}\" ]]; then echo $idx; fi; done)]}}; "
        "md5checksum=$(md5sum \"{input.thing}\" | awk '{{print $1}}'); "
        "echo \"${{id}};{input.thing};Created output;$(stat -c '%.19z' {input.thing});${{md5checksum}};MeBaPiNa;"
        "K-mer: taxonomic classification.\" "
        ">> \"${{report_log}}\"; "
        "touch {output.dummy}"

rule report_retax_kmermap:
    input:
        thing="{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/{reftype}.kreport2"
    output:
        dummy=temp("{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_retax_kmermap.report")
    params:
        "report_log=" + LOGS + "; ",
        "all_IDs=( $(echo \"" + " ".join(SAMPLES.values()) + "\") ); ",
        "all_barcs=( $(echo \"" + " ".join(SAMPLES.keys()) + "\") ); "
    shell:
        "{params}"
        "id=${{all_IDs[$(for idx in \"${{!all_barcs[@]}}\"; do "
        "if [[ \"${{all_barcs[$idx]}}\" = \"{wildcards.barc}\" ]]; then echo $idx; fi; done)]}}; "
        "md5checksum=$(md5sum \"{input.thing}\" | awk '{{print $1}}'); "
        "echo \"${{id}};{input.thing};Created output;$(stat -c '%.19z' {input.thing});${{md5checksum}};MeBaPiNa;"
        "K-mer: abundance reestimation.\" "
        ">> \"${{report_log}}\"; "
        "touch {output.dummy}"

rule report_counttax_kmermap:
    input:
        thing="{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/kmer.counttaxlist"
    output:
        dummy=temp("{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_counttax_kmermap.report")
    params:
        "report_log=" + LOGS + "; ",
        "all_IDs=( $(echo \"" + " ".join(SAMPLES.values()) + "\") ); ",
        "all_barcs=( $(echo \"" + " ".join(SAMPLES.keys()) + "\") ); "
    shell:
        "{params}"
        "id=${{all_IDs[$(for idx in \"${{!all_barcs[@]}}\"; do "
        "if [[ \"${{all_barcs[$idx]}}\" = \"{wildcards.barc}\" ]]; then echo $idx; fi; done)]}}; "
        "md5checksum=$(md5sum \"{input.thing}\" | awk '{{print $1}}'); "
        "echo \"${{id}};{input.thing};Created output;$(stat -c '%.19z' {input.thing});${{md5checksum}};MeBaPiNa;"
        "K-mer: file conversion.\" "
        ">> \"${{report_log}}\"; "
        "touch {output.dummy}"


## CALIBRATION STRAIN ##
########################
