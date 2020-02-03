########################
## SEQUENCING SUMMARY ##
########################

rule sorting_seqsum_barc:
    input:
        "01_processed_data/01_basecalling/{run}/sequencing_summary.txt"
    output:
        "01_processed_data/01_basecalling/{run}/sequencing_summary/sequencing_summary_sorted.txt"
    log:
        "01_processed_data/01_basecalling/{run}/sequencing_summary/MeBaPiNa_sorting.log"
    benchmark:
        "01_processed_data/01_basecalling/{run}/sequencing_summary/MeBaPiNa_sorting.benchmark.tsv"
    threads:
        4
    shell:
        "barc_col=$( awk '{{ header=$0; " ## save header string and ...
        "for(i;i<=NF;i++){{if($i==\"barcode_arrangement\"){{ barc_col=i }}}} }}; " ## ...find column with barcode name
        "{{ print header > \"{output}\"; print barc_col; exit 0 }}' {input} ) > {log} 2>&1; " ## write header to file and print barcode column number to save to variable, exit after first line
        "tail -n +2 {input} | "
        "sort -V -S1G --parallel={threads} -k$barc_col,$barc_col -k7,7 " ## sort by barcode column and start time (Note, start time column might have another index in later versions)
        ">> {output} 2>> {log}"

rule splitting_seqsum_barc:
    input:
        "01_processed_data/01_basecalling/{run}/sequencing_summary/sequencing_summary_sorted.txt"
    output:
        directory("01_processed_data/01_basecalling/{run}/sequencing_summary/split")
    log:
        "01_processed_data/01_basecalling/{run}/sequencing_summary/MeBaPiNa_splitting.log"
    benchmark:
        "01_processed_data/01_basecalling/{run}/sequencing_summary/MeBaPiNa_splitting.benchmark.tsv"
    shell:
        "mkdir -p {output} > {log} 2>&1; "
        "awk 'BEGIN{{cali_col=0}}; " ## falsify calibration strain column detection
        "NR==1{{ header=$0; " ## save header string (first line) and...
        "for(i;i<=NF;i++){{ if($i==\"barcode_arrangement\"){{ barc_col=i }}; " ## ...find column with barcode name and...
        "if($i==\"calibration_strand_genome\"){{ cali_col=i; " ## ...detect the presence of calibration strands and...
        "print header > \"{output}\" \"/lambda.txt\" }} }} }}; " ## ...when detected, write header to calibration specific file
        "NR!=1{{ " ## for all other lines
        "if(firstencounter[$barc_col]==0){{ " ## when the barcode is encountered the first time (no file exsist, yet)...
        "print header > \"{output}\" \"/\" $barc_col \".txt\"; " ## ...write the header string to the barcode specific file
        "firstencounter[$barc_col]++ }}; " ## set the firstencounter for this barcode to false (!=0)
        "if( cali_col != 0 && $barc_col == \"unclassified\" && $cali_col ~ /Lambda_3.6kb/){{ print $0 > \"{output}\" \"/lambda.txt\" }}"
        "else{{ print $0 > \"{output}\" \"/\" $barc_col \".txt\" }} " ## print the current line to the corresponding barcode specific file
        "}}' {input} >> {log} 2>&1"

rule downsampling_seqsum: #!#
    input:
        "01_processed_data/01_basecalling/{run}/sequencing_summary/sequencing_summary_sorted.txt"
    output:
        "01_processed_data/01_basecalling/{run}/sequencing_summary/sequencing_summary_downsampled.txt"
    log:
        "01_processed_data/01_basecalling/{run}/sequencing_summary/MeBaPiNa_downsampling.log"
    benchmark:
        "01_processed_data/01_basecalling/{run}/sequencing_summary/MeBaPiNa_downsampling.benchmark.tsv"
    shell:
        "cat {input} | (read -r; printf \"%s\\n\" \"$REPLY\"; "
        "awk -v seed=$RANDOM 'BEGIN{{ prnt=-4; nr=100; srand(seed) }}; " ## "falsify" print flag, set nr for fraction and set random seed. Note: downsampling is a fraction here not a number of reads
        "{{ rhundr=1+int(rand()*nr) }}; " ## get a random number between 1 and nr
        "rhundr==nr') " ## if the number is nr (by chance of 1/nr) print the line
        ">> {output} 2> {log}"
