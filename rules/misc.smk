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

################
## REFERENCES ##
################

## SILVA ##
###########

rule download_reffiles:
    output:
        fasta="METADATA/Reference_Sequences/silva/reference.fasta",
        ncbimap="METADATA/Reference_Sequences/silva/ncbimap.txt",
        slvmap="METADATA/Reference_Sequences/silva/slvmap.txt",
        taxlist="METADATA/Reference_Sequences/silva/taxlist.txt"
        # fasta2="{reference}/kraken2/species_tmp/library/silva.fna",
        # fasta3="{reference}/kraken2/genus_tmp/library/silva.fna"
    log:
        "METADATA/Reference_Sequences/silva/MeBaPiNa_reffiles.log"
    benchmark:
        "METADATA/Reference_Sequences/silva/MeBaPiNa_reffiles.benchmark.tsv"
    shell:
        "bash Pineline/MeBaPiNa/scripts/download_silva.sh METADATA/Reference_Sequences/silva >> {log} 2>&1"

rule krona_reffile:
    output:
        directory("METADATA/Reference_Sequences/ncbi/krona")
    log:
        "METADATA/Reference_Sequences/ncbi/MeBaPiNa_reffiles.log"
    benchmark:
        "METADATA/Reference_Sequences/ncbi/MeBaPiNa_reffiles.benchmark.tsv"
    conda:
        "../envs/krona.yml"
    shell:
        "ktUpdateTaxonomy.sh {output} >> {log} 2>&1"

rule construct_reffiles:
    input:
        ncbimap="METADATA/Reference_Sequences/silva/ncbimap.txt",
        slvmap="METADATA/Reference_Sequences/silva/slvmap.txt",
        taxlist="METADATA/Reference_Sequences/silva/taxlist.txt",
        ncbikrona="METADATA/Reference_Sequences/ncbi/krona"
    output:
        krak_S=directory("METADATA/Reference_Sequences/silva/kraken2/species_tmp"),
        krak_G=directory("METADATA/Reference_Sequences/silva/kraken2/genus_tmp"),
        krona_S=directory("METADATA/Reference_Sequences/silva/krona/species"),
        krona_G=directory("METADATA/Reference_Sequences/silva/krona/genus"),
        taxlist_S="METADATA/Reference_Sequences/silva/krona/species/taxlist.txt",
        taxlist_G="METADATA/Reference_Sequences/silva/krona/genus/taxlist.txt"
    log:
        "METADATA/Reference_Sequences/silva/MeBaPiNa_reffiles.log"
    benchmark:
        "METADATA/Reference_Sequences/silva/MeBaPiNa_reffiles.benchmark.tsv"
    script:
        "../scripts/construct_reffiles.py"

## KRAKEN2 DATABASE ##
######################

rule building_database_fromreffiles:
    input:
        "METADATA/Reference_Sequences/{reference}/kraken2/{reftype}_tmp"
    output:
        directory("METADATA/Reference_Sequences/{reference}/kraken2/{reftype}")
    log:
        "METADATA/Reference_Sequences/{reference}/kraken2/MeBaPiNa_{reftype}.log"
    benchmark:
        "METADATA/Reference_Sequences/{reference}/kraken2/MeBaPiNa_{reftype}.benchmark.tsv"
    conda:
        "../envs/kraken2.yml"
    threads:
        8
    params:
        "35" ## k-mer length
    shell:
        "mv {input} {output}; "
        # "kraken2-build --threads {threads} --download-taxonomy --skip-maps --db {output} > {log} 2>&1; "
        # "kraken2-build --threads {threads} --download-library bacteria --no-masking --db {output} >> {log} 2>&1; "
        # # "kraken2-build --threads {threads} --download-library archaea --no-masking --db {output} >> {log} 2>&1; "
        "kraken2-build --threads {threads} --kmer-len {params} --build "
        "--db {output} "
        ">> {log} 2>&1; "
        # "kraken2-build --clean --db {output} >> {log} 2>&1"
        # "kraken2-build --threads {threads} " ## unfortunately --kmer-len {params} is ignored by when --special is used
        # "--special {wildcards.reference} --db {output} > {log}; " ## reference can be one of "greengenes", "silva", "rdp"
        "bracken-build -t {threads} -k {params} -l 1451 " ## 1451 ismedian read length after filtering in 20191007_1559_MN31344_FAK76605_2bf006ff
        "-d {output} "
        ">> {log} 2>&1; "
        "kraken2-build --clean --db {output} >> {log} 2>&1"

rule building_database_alone:
    output:
        directory("METADATA/Reference_Sequences/{reference}/kraken2/{reftype}") ## reftype has no incluence on this rule, but is used for silva species database (rule building_database_fromreffiles)
    log:
        "METADATA/Reference_Sequences/{reference}/kraken2/MeBaPiNa_{reftype}.log"
    benchmark:
        "METADATA/Reference_Sequences/{reference}/kraken2/MeBaPiNa_{reftype}.benchmark.tsv"
    conda:
        "../envs/kraken2.yml"
    threads:
        8
    params:
        "35" ## k-mer length
    shell:
        "kraken2-build --threads {threads} " ## unfortunately --kmer-len {params} is ignored by when --special is used
        "--special {wildcards.reference} --db {output} " ## reference can be one of "greengenes", "silva", "rdp"
        "> {log}; "
        "bracken-build -t {threads} -k {params} -l 1451 -d {output} " ## 1451 ismedian read length after filtering in 20191007_1559_MN31344_FAK76605_2bf006ff
        ">> {log}; "
        "rm -rf {output}/data "
        "kraken2-build --clean --db {output} >> {log} 2>&1"

ruleorder: building_database_fromreffiles > building_database_alone

## MINIMA2 INDEX ##
###################

rule indexing_reference:
    input:
        "METADATA/Reference_Sequences/{reference}/reference.fasta"
    output:
        "METADATA/Reference_Sequences/{reference}/reference.mmi"
    log:
        "METADATA/Reference_Sequences/{reference}/MeBaPiNa_indexing.log"
    benchmark:
        "METADATA/Reference_Sequences/{reference}/MeBaPiNa_indexing.benchmark.tsv"
    params:
        "-x map-ont" ## naopore specific
    conda:
        "../envs/minimap2.yml"
    threads:
        2
    shell:
        "minimap2 -t {threads} {params} -d {output} {input} > {log} 2>&1"
