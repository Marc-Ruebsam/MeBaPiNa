#####################
## FILE CONVERSION ##
#####################

rule sam_to_bam:
    input:
        "01_processeddata/{run}/{align}/{barc}_alignment.sam"
    output:
        bam="01_processeddata/{run}/{align}/{barc}_alignment_sorted.bam",
        bai="01_processeddata/{run}/{align}/{barc}_alignment_sorted.bam.bai"
    log:
        "01_processeddata/{run}/{align}/{barc}_MeBaPiNa_sam_bam.log"
    benchmark:
        "01_processeddata/{run}/{align}/{barc}_MeBaPiNa_sam_bam.benchmark.tsv"
    conda:
        "../envs/samtools.yml"
    threads:
        2
    shell:
        "samtools sort --threads {threads} -o {output.bam} {input} > {log} 2>&1; "
        "samtools index -@ {threads} {output.bam} >> {log} 2>&1"

########################
## SEQUENCING SUMMARY ##
########################

rule find_reads_in_fastq:
    input:
        "01_processeddata/{run}/basecall/pass"
    output:
        "01_processeddata/{run}/basecall/find_reads_in_fastq.txt"
    log:
        "01_processeddata/{run}/basecall/MeBaPiNa_find_fastq_in_fast5.log"
    benchmark:
        "01_processeddata/{run}/basecall/MeBaPiNa_find_fastq_in_fast5.benchmark.tsv"
    conda:
        "../envs/find_fastq_in_fast5.yml"
    script:
        "../scripts/find_fastq_in_fast5.py"

rule split_seqsum_barc:
    input:
        "01_processeddata/{run}/basecall/sequencing_summary/sequencing_summary_sorted.txt"
    output:
        directory("01_processeddata/{run}/basecall/sequencing_summary/barcode")
    log:
        "01_processeddata/{run}/basecall/sequencing_summary/MeBaPiNa_split_seqsum_barc.log"
    benchmark:
        "01_processeddata/{run}/basecall/sequencing_summary/MeBaPiNa_split_seqsum_barc.benchmark.tsv"
    shell:
        "mkdir -p {output} > {log} 2>&1; "
        "awk 'BEGIN{{cali_col=0}}; " ## falsify calibration strain column detection
        "NR==1{{ header=$0; " ## save header string (first line) and...
        "for(i;i<=NF;i++){{ if($i==\"barcode_arrangement\"){{ barc_col=i }}; " ## ...find column with barcode name and...
        "if($i==\"calibration_strand_genome\"){{ cali_col=i; " ## ...detect the presence of calibration strands and...
        "print header > \"{output}\" \"/sequencing_summary_lambda.txt\" }} }} }}; " ## ... write header to calibration specific file
        "NR!=1{{ " ## for all other lines
        "if(firstencounter[$barc_col]==0){{ " ## when the barcode is encountered the first time (no file exsist, yet)...
        "print header > \"{output}\" \"/sequencing_summary_\" $barc_col \".txt\"; " ## ...write the header string to the barcode specific file
        "firstencounter[$barc_col]++ }}; " ## set the firstencounter for this string to false (!=0)
        "if( cali_col != 0 && $barc_col == \"unclassified\" && $cali_col ~ /Lambda_3.6kb/){{ print $0 > \"{output}\" \"/sequencing_summary_lambda.txt\" }}"
        "else{{ print $0 > \"{output}\" \"/sequencing_summary_\" $barc_col \".txt\" }} " ## print the current line to the corresponding barcode specific file
        "}}' {input} >> {log} 2>&1"

rule sort_seqsum:
    input:
        "01_processeddata/{run}/basecall/sequencing_summary.txt"
    output:
        "01_processeddata/{run}/basecall/sequencing_summary/sequencing_summary_sorted.txt"
    log:
        "01_processeddata/{run}/basecall/sequencing_summary/MeBaPiNa_sort_seqsum_barc.log"
    benchmark:
        "01_processeddata/{run}/basecall/sequencing_summary/MeBaPiNa_sort_seqsum_barc.benchmark.tsv"
    threads:
        4
    shell:
        "barc_col=$( awk '{{ header=$0; " ## save header string and ...
        "for(i;i<=NF;i++){{if($i==\"barcode_arrangement\"){{ barc_col=i }}}} }}; " ## ...find column with barcode name
        "{{ print header > \"{output}\"; print barc_col; exit 0 }}' {input} ) > {log} 2>&1; " ## write header to file and print barcode column number to save to variable, exit after first line
        "tail -n +2 {input} | "
        "sort -V -S1G --parallel={threads} -k$barc_col,$barc_col -k7,7 " ## sort by barcode column and start time (Note, start time column might have another index in later versions)
        ">> {output} 2>> {log}"

rule downsample_seqsum:
    input:
        "01_processeddata/{run}/basecall/sequencing_summary/sequencing_summary_sorted.txt"
    output:
        "01_processeddata/{run}/basecall/sequencing_summary/sequencing_summary_downsampled.txt"
    log:
        "01_processeddata/{run}/basecall/sequencing_summary/MeBaPiNa_downsample_seqsum.log"
    benchmark:
        "01_processeddata/{run}/basecall/sequencing_summary/MeBaPiNa_downsample_seqsum.benchmark.tsv"
    shell:
        "cat {input} | (read -r; printf \"%s\\n\" \"$REPLY\"; "
        "awk -v seed=$RANDOM 'BEGIN{{ prnt=-4; nr=100; srand(seed) }}; " ## "falsify" print flag, set nr for fraction and set random seed. Note: downsampling is a fraction here not a number of reads
        "{{ rhundr=1+int(rand()*nr) }}; " ## get a random number between 1 and nr
        "rhundr==nr') " ## if the number is nr (by chance of 1/nr) print the line
        ">> {output} 2> {log}"

############################
## CHECKPOINT AGGREGATION ##
############################

def input_aggregate(wildcards):
    from os import listdir
    ## "pass" directory
    basecall_dir = checkpoints.guppy.get(run=wildcards.run).output[0]
    ## directory names within "pass" directory
    barcs = listdir(basecall_dir)
    ## retain only strings containing "barcode"
    barcs = [i for i in barcs if "barcode" in i]
    ## create file names with barcodes
    barc_input = expand("02_analysis/{run}/align/{barc}_pycoqc/pycoQC_report.json",
        run=wildcards.run,
        barc=barcs)
    barc_input.sort()
    return barc_input

rule aggr_align_barc:
    input:
        input_aggregate ## if it keeps returning only the pass folder, some but not all of the output and input of the guppy rule exists. Create it, even if its only an empty folder to fix the problem.
    output:
        "02_analysis/{run}/align/MeBaPiNa_barcode_aggregation.txt"
    shell:
        "echo {input} > {output}"
