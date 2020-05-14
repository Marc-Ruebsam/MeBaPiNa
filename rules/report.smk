##############
## PLOTTING ##
##############

## BASECALLING ##
#################

rule report_moving_raw:
    input:
        "{tmp}00_raw_data/{run}/fast5"
    output:
        temp("{tmp}00_raw_data/{run}/MeBaPiNa_moving_raw.report")
    params:
        " ".join(SAMPLES.values())
    shell:
        "all_IDs=( $(echo \"{params}\") ); "
        "for id in ${{all_IDs[@]}}; do echo "
        ## "Sample name;File/directory;Completion date;Checksum;Performed by;Description"
        "\"${{id}};{input};$(stat -c %y {input});NA;MeBaPiNa;Moved raw fast5 files to path ready for analysis.\" "
        ">> {output}; done"

rule report_basecalling_raw:
    input:
        list(filter(None,[
        "{tmp}01_processed_data/01_basecalling/{run}/sequencing_summary.txt",
        directory("{tmp}01_processed_data/01_basecalling/{run}/pass"),
        (directory("{tmp}01_processed_data/01_basecalling/{run}/calibration_strands") if LAM_DCS else "") ## evaluation of Lambda calibration strands only when specified
        ]))
    output:
        temp("{tmp}00_raw_data/{run}/MeBaPiNa_basecalling_raw.report")
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

rule report_trimming_basecalled:
    input:
        "{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/trimmed.fastq"
    output:
        temp("{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/MeBaPiNa_trimming_basecalled.report")
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

rule report_filtering_trimmed:
    input:
        "{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/filtered.fastq"
    output:
        temp("{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/MeBaPiNa_filtering_trimmed.report")
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

## ALIGN ##
###########

## K-MER MAPPING ##
###################

## CALIBRATION STRAIN ##
########################
