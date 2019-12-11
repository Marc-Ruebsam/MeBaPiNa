##############################
## TRIMM FILTER DEMULTIPLEX ##
##############################

rule nanofilt:
    input:
        fastq="01_processeddata/{run}/basecall/pass/{barc}",
        seqsum="01_processeddata/{run}/basecall/sequencing_summary.txt"
    output:
        "01_processeddata/{run}/filter/{barc}.fastq.gz"
    log:
        "01_processeddata/{run}/filter/{barc}_MeBaPiNa_nanofilt.log"
    benchmark:
        "01_processeddata/{run}/filter/{barc}_MeBaPiNa_nanofilt.benchmark.tsv"
    conda:
        "../envs/nanopack.yml"
    params:
        ("--length " + config["guppy"]["len_min"]),
        ("--maxlength " + config["guppy"]["len_max"]),
        ("--quality " + config["guppy"]["q_min"])
    shell:
        "gunzip -c $(find {input.fastq} -type f -name \"*.fastq.gz\") | "
        "NanoFilt {params} --summary {input.seqsum} "
        "--logfile {log} "
        "| gzip > {output}"

rule porechop:
    input:
        "01_processeddata/{run}/filter/{barc}.fastq.gz"
    output:
        "01_processeddata/{run}/trim/{barc}.fastq.gz"
    log:
        "01_processeddata/{run}/trim/{barc}_MeBaPiNa_porechop.log"
    benchmark:
        "01_processeddata/{run}/trim/{barc}_MeBaPiNa_porechop.benchmark.tsv"
    conda:
        "../envs/porechop.yml"
    threads:
        config["machine"]["cpu"]
    params:
        "--format fastq",
        "--check_reads 1000", ## This many reads will be aligned to all possible adapters to determine which adapter sets are present (default: 10000)
        "--barcode_threshold 75.0", ## A read must have at least this percent identity to a barcode to be binned (default: 75.0)
        "--barcode_diff 5.0", ## If the difference between a read's best barcode identity and its second-best barcode identity is less than this value, it will not be put in a barcode bin (to exclude cases which are too close to call) (default: 5.0)
        # "--require_two_barcodes", ## Reads will only be put in barcode bins if they have a strong match for the barcode on both their start and end (default: a read can be binned with a match at its start or end)
        "--adapter_threshold 90.0", ## An adapter set has to have at least this percent identity to be labelled as present and trimmed off (0 to 100) (default: 90.0)
        "--end_threshold 75.0", ## Adapters at the ends of reads must have at least this percent identity to be removed (0 to 100) (default: 75.0)
        "--extra_end_trim 2", ## This many additional bases will be removed next to adapters found at the ends of reads (default: 2)
        "--verbosity 1" ##  Level of progress information: 0 = none, 1 = some, 2 = lots, 3 = full - output will go to stdout if reads are saved to a file and stderr if reads are printed to stdout (default: 1)
    shell:
        "mkdir -p 01_processeddata/{wildcards.run}/trim/{wildcards.barc} > {log} 2>&1; "
        "porechop --threads {threads} {params} --input {input} --barcode_dir 01_processeddata/{wildcards.run}/trim/{wildcards.barc} > {log} 2>&1; "
        "barcode_string={wildcards.barc} > {log} 2>&1; "
        "gzip --stdout 01_processeddata/{wildcards.run}/trim/{wildcards.barc}/BC${{barcode_string: -2}}.fastq > 01_processeddata/{wildcards.run}/trim/{wildcards.barc}.fastq.gz 2> {log}"

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
        2
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
        "01_processeddata/{run}/basecall/sequencing_summary/sequencing_summary_sorted.txt"
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
        "{{ print header > \"{output}\"; print barc_col; exit 0 }}' ); " ## write header to file and print barcode column number to save to variable, exit after first line
        "sort -V --parallel={threads} -k$barc_col,$barc_col -k7,7 <( tail -n +2 {input} )" ## sort by barcode column and start time (Note, start time column might have another index in later versions)
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
    ## "pass" directory
    basecall_dir = checkpoints.guppy.get(run=wildcards.run).output[0]
    ## directory names within "pass" directory
    barcs = listdir(basecall_dir)
    ## retain only strings containing "barcode"
    barcs = [i for i in barcs if "barcode" in i]
    ## create file names with barcodes
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
