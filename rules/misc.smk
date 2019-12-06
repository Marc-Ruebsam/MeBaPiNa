########################################
## FILE CONVERSION AND CONCATENATION ##
########################################

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
    version:
        subprocess.check_output("samtools --version | awk 'NR==1{print $NF}'", shell=True)
    conda:
        "../envs/samtools.yml"
    threads:
        config["machine"]["cpu"]
    shell:
        "samtools sort --threads {threads} -o {output.bam} {input} > {log} 2>&1; "
        "samtools index -@ {threads} {output.bam} > {log} 2>&1"

rule fasq_pipe:
    input:
        "01_processeddata/{run}/basecall/pass"
    output:
        temp("02_analysis/{run}/basecall/nanoqc/pipe.fastq") ## pipe didn't work neighter did named pipes (no fastq file extension)
    log:
        "02_analysis/{run}/basecall/nanoqc/MeBaPiNa_fasq_pipe.log"
    benchmark:
        "02_analysis/{run}/basecall/nanoqc/MeBaPiNa_fasq_pipe.benchmark.tsv"
    shell:
        "zcat $(find {input} -type f -name \"*.fastq.gz\") |" ## concatenate all fastq files
        "awk -v seed=$RANDOM 'BEGIN{{ prnt=-4; nr=100; srand(seed) }}; " ## "falsify" print flag, set nr for fraction and set random seed. Note: downsampling is a fraction here not a number of reads
        "NR%4==1{{ " ## for every fourth line -> all headers
        "rhundr=1+int(rand()*nr); " ## get a random number between 1 and nr
        "if(rhundr==nr)" ## if the number is nr (by chance of 1/nr)
        "{{ prnt=NR }} " ## set prnt flag to current line number to... 
        "}}; "
        "NR<prnt+4' " ## ...print this and the next three lines
        ">> {output} 2> {log}"

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
        "01_processeddata/{run}/basecall/sequencing_summary.txt"
    output:
        directory("01_processeddata/{run}/basecall/sequencing_summary/barcode")
    log:
        "01_processeddata/{run}/basecall/sequencing_summary/MeBaPiNa_split_seqsum_barc.log"
    benchmark:
        "01_processeddata/{run}/basecall/sequencing_summary/MeBaPiNa_split_seqsum_barc.benchmark.tsv"
    shell:
        "mkdir {output}; "
        "awk 'NR==1{{ header=$0; " ## save header string (first line) and ...
        "for(i;i<=NF;i++){{if($i==\"barcode_arrangement\"){{ barc_col=i }}}} }}; " ## ...find column with barcode name
        "NR!=1{{ " ## for all other lines
        "if(firstencounter[$barc_col]==0){{ " ## when the barcode is encountered the first time (no file exsist, yet)...
        "print header > \"{output}\" \"/sequencing_summary_\" $barc_col \".txt\"; " ## ...write the header string to the barcode specific file
        "firstencounter[$barc_col]++ }}; " ## set the firstencounter for this string to false (!=0)
        "print $0 > \"{output}\" \"/sequencing_summary_\" $barc_col \".txt\" " ## print the current line to the corresponding barcode specific file
        "}}' {input} > {log} 2>&1"

rule sort_seqsum_barc:
    input:
        "01_processeddata/{run}/basecall/sequencing_summary.txt"
    output:
        "01_processeddata/{run}/basecall/sequencing_summary/sequencing_summary_sorted.txt"
    log:
        "01_processeddata/{run}/basecall/sequencing_summary/MeBaPiNa_sort_seqsum_barc.log"
    benchmark:
        "01_processeddata/{run}/basecall/sequencing_summary/MeBaPiNa_sort_seqsum_barc.benchmark.tsv"
    threads:
        config["machine"]["cpu"]
    shell:
        "cat {input} | (read -r; printf \"%s\\n\" \"$REPLY\"; sort -V --parallel={threads} -k"
        "$(awk '{{ for(i;i<=NF;i++){{if($i==\"barcode_arrangement\"){{ print i }}}}; exit 0 }}' {input})," ## finds the "barcode_arrangement" column...
        "$(awk '{{ for(i;i<=NF;i++){{if($i==\"barcode_arrangement\"){{ print i }}}}; exit 0 }}' {input}) -k7,7) " ## ...twice
        ">> {output} 2> {log}"

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
    basecall_dir = checkpoints.guppy.get(run=wildcards.run).output[0]
    barcs = listdir(basecall_dir)
    barc_dirs = expand("02_analysis/{run}/align/{barc}_pycoqc/pycoQC_report.json",
        run=wildcards.run,
        barc=barcs)
    return barc_dirs

rule aggr_align_barc:
    input:
        input_aggregate ## if it keeps returning only the pass folder, some but not all of the output and input of the guppy rule exists. Create it, even if its only an empty folder to fix the problem.
    output:
        "02_analysis/{run}/align/MeBaPiNa_barcode_aggregation.txt"
    shell:
        "echo {input} > {output}"
