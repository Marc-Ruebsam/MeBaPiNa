########################
## SEQUENCING SUMMARY ##
########################

## SORTING ##
#############

rule sorting_seqsum_barc:
    input:
        "{tmp}01_processed_data/01_basecalling/{run}/sequencing_summary.txt"
    output:
        "{tmp}01_processed_data/01_basecalling/{run}/sequencing_summary/sequencing_summary_sorted.txt"
    log:
        "{tmp}01_processed_data/01_basecalling/{run}/sequencing_summary/MeBaPiNa_sorting.log"
    benchmark:
        "{tmp}01_processed_data/01_basecalling/{run}/sequencing_summary/MeBaPiNa_sorting.benchmark.tsv"
    threads:
        4
    shell:
        "barc_col=$( awk '{{ header=$0; " ## save header string and ...
        "for(i;i<=NF;i++){{if($i==\"barcode_arrangement\"){{ barc_col=i }}}} }}; " ## ...find column with barcode name
        "{{ print header > \"{output}\"; print barc_col; exit 0 }}' {input} ) > {log} 2>&1; " ## write header to file and print barcode column number to save to variable, exit after first line
        "tail -n +2 {input} | "
        "sort -V -S1G --parallel={threads} -k$barc_col,$barc_col -k7,7 " ## sort by barcode column and start time (Note, start time column might have another index in later versions)
        ">> {output} 2>> {log}"

## SPLITTING ##
###############

rule splitting_seqsum_barc:
    input:
        "{tmp}01_processed_data/01_basecalling/{run}/sequencing_summary/sequencing_summary_sorted.txt"
    output:
        directory("{tmp}01_processed_data/01_basecalling/{run}/sequencing_summary/split")
    log:
        "{tmp}01_processed_data/01_basecalling/{run}/sequencing_summary/MeBaPiNa_splitting.log"
    benchmark:
        "{tmp}01_processed_data/01_basecalling/{run}/sequencing_summary/MeBaPiNa_splitting.benchmark.tsv"
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

## DOWNSAMPLING ##
##################

rule downsampling_seqsum: #!#
    input:
        "{tmp}01_processed_data/01_basecalling/{run}/sequencing_summary/sequencing_summary_sorted.txt"
    output:
        "{tmp}01_processed_data/01_basecalling/{run}/sequencing_summary/sequencing_summary_downsampled.txt"
    log:
        "{tmp}01_processed_data/01_basecalling/{run}/sequencing_summary/MeBaPiNa_downsampling.log"
    benchmark:
        "{tmp}01_processed_data/01_basecalling/{run}/sequencing_summary/MeBaPiNa_downsampling.benchmark.tsv"
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
        fasta="{tmp}METADATA/Reference_Sequences/silva/reference.fasta",
        fasta_align="{tmp}METADATA/Reference_Sequences/silva/reference_aligned.fasta",
        ncbimap="{tmp}METADATA/Reference_Sequences/silva/ncbimap.txt",
        slvmap="{tmp}METADATA/Reference_Sequences/silva/slvmap.txt",
        taxlist="{tmp}METADATA/Reference_Sequences/silva/taxlist.txt",
        taxtre="{tmp}METADATA/Reference_Sequences/silva/taxlist.tre"
    log:
        "{tmp}METADATA/Reference_Sequences/silva/MeBaPiNa_reffiles.log"
    benchmark:
        "{tmp}METADATA/Reference_Sequences/silva/MeBaPiNa_reffiles.benchmark.tsv"
    shell:
        "bash {wildcards.tmp}Pineline/MeBaPiNa/scripts/download_silva.sh {wildcards.tmp}METADATA/Reference_Sequences/silva > {log} 2>&1"

rule krona_reffile:
    output:
        "{tmp}METADATA/Reference_Sequences/ncbi/krona/taxonomy.tab"
    log:
        "{tmp}METADATA/Reference_Sequences/ncbi/MeBaPiNa_reffiles.log"
    benchmark:
        "{tmp}METADATA/Reference_Sequences/ncbi/MeBaPiNa_reffiles.benchmark.tsv"
    conda:
        "../envs/krona.yml"
    shell:
        "ktUpdateTaxonomy.sh {output} > {log} 2>&1"

rule construct_reffiles:
    input:
        taxlist="{tmp}METADATA/Reference_Sequences/silva/taxlist.txt",
        slvmap="{tmp}METADATA/Reference_Sequences/silva/slvmap.txt",
        ncbimap="{tmp}METADATA/Reference_Sequences/silva/ncbimap.txt",
        ncbikrona="{tmp}METADATA/Reference_Sequences/ncbi/krona/taxonomy.tab"
    output:
        kraknames_S="{tmp}METADATA/Reference_Sequences/silva/kraken2/species/taxonomy/names.dmp",
        kraknodes_S="{tmp}METADATA/Reference_Sequences/silva/kraken2/species/taxonomy/nodes.dmp",
        krakseq2tax_S="{tmp}METADATA/Reference_Sequences/silva/kraken2/species/seqid2taxid.map",
        kraknames_G="{tmp}METADATA/Reference_Sequences/silva/kraken2/genus/taxonomy/names.dmp",
        kraknodes_G="{tmp}METADATA/Reference_Sequences/silva/kraken2/genus/taxonomy/nodes.dmp",
        krakseq2tax_G="{tmp}METADATA/Reference_Sequences/silva/kraken2/genus/seqid2taxid.map",
        kronataxtab_S="{tmp}METADATA/Reference_Sequences/silva/krona/species/taxonomy.tab",
        kronataxlist_S="{tmp}METADATA/Reference_Sequences/silva/krona/species/taxlist.txt",
        kronaseq2tax_S="{tmp}METADATA/Reference_Sequences/silva/krona/species/seqid2taxid.map",
        kronataxtab_G="{tmp}METADATA/Reference_Sequences/silva/krona/genus/taxonomy.tab",
        kronataxlist_G="{tmp}METADATA/Reference_Sequences/silva/krona/genus/taxlist.txt",
        kronaseq2tax_G="{tmp}METADATA/Reference_Sequences/silva/krona/genus/seqid2taxid.map"
    log:
        "{tmp}METADATA/Reference_Sequences/silva/MeBaPiNa_reffiles.log"
    benchmark:
        "{tmp}METADATA/Reference_Sequences/silva/MeBaPiNa_reffiles.benchmark.tsv"
    conda:
        "../envs/python.yml"
    script:
        "../scripts/construct_reffiles.py"

## KRAKEN2 DATABASE ##
######################

rule building_database_fromreffiles:
    input:
        fastaT="{tmp}METADATA/Reference_Sequences/{reference}/reference_thymine.fasta",
        kraknames="{tmp}METADATA/Reference_Sequences/{reference}/kraken2/{reftype}/taxonomy/names.dmp",
        kraknodes="{tmp}METADATA/Reference_Sequences/{reference}/kraken2/{reftype}/taxonomy/nodes.dmp",
        krakseq2tax="{tmp}METADATA/Reference_Sequences/{reference}/kraken2/{reftype}/seqid2taxid.map"
    output:
        krakdb="{tmp}METADATA/Reference_Sequences/{reference}/kraken2/{reftype}/database.kraken",
        krakhash="{tmp}METADATA/Reference_Sequences/{reference}/kraken2/{reftype}/hash.k2d",
        krakopts="{tmp}METADATA/Reference_Sequences/{reference}/kraken2/{reftype}/opts.k2d",
        kraktaxo="{tmp}METADATA/Reference_Sequences/{reference}/kraken2/{reftype}/taxo.k2d",
        brakdb="{tmp}METADATA/Reference_Sequences/{reference}/kraken2/{reftype}/database1451mers.kraken",
        brakdist="{tmp}METADATA/Reference_Sequences/{reference}/kraken2/{reftype}/database1451mers.kmer_distrib"
    log:
        "{tmp}METADATA/Reference_Sequences/{reference}/kraken2/MeBaPiNa_{reftype}.log"
    benchmark:
        "{tmp}METADATA/Reference_Sequences/{reference}/kraken2/MeBaPiNa_{reftype}.benchmark.tsv"
    conda:
        "../envs/kraken2.yml"
    threads:
        6
    params:
        "35" ## k-mer length
    shell:
        "out_dir={output.krakdb}; out_dir=\"${{out_dir/database.kraken/}}\" > {log} 2>&1; "
        "mkdir \"${{out_dir}}library\" >> {log} 2>&1; "
        "cp \"{input.fastaT}\" \"${{out_dir}}library/library.fna\" >> {log} 2>&1; "
        # "kraken2-build --threads {threads} --download-taxonomy --skip-maps --db {output} > {log} 2>&1; "
        # "kraken2-build --threads {threads} --download-library bacteria --no-masking --db {output} >> {log} 2>&1; "
        # # "kraken2-build --threads {threads} --download-library archaea --no-masking --db {output} >> {log} 2>&1; "
        "kraken2-build --threads {threads} --kmer-len {params} --build "
        "--db ${{out_dir}} "
        ">> {log} 2>&1; "
        # "kraken2-build --clean --db {output} >> {log} 2>&1"
        # "kraken2-build --threads {threads} " ## unfortunately --kmer-len {params} is ignored by when --special is used
        # "--special {wildcards.reference} --db {output} > {log}; " ## reference can be one of "greengenes", "silva", "rdp"
        "bracken-build -t {threads} -k {params} -l 1451 " ## 1451 ismedian read length after filtering in 20191007_1559_MN31344_FAK76605_2bf006ff
        "-d ${{out_dir}} "
        ">> {log} 2>&1; "
        "kraken2-build --clean --db ${{out_dir}} >> {log} 2>&1"

rule building_database_alone:
    output:
        krakdb="{tmp}METADATA/Reference_Sequences/{reference}/kraken2/{reftype}/database.kraken", ## reftype wildcard has no influence on this rule, but is used for silva species database (rule building_database_fromreffiles)
        krakhash="{tmp}METADATA/Reference_Sequences/{reference}/kraken2/{reftype}/hash.k2d",
        krakopts="{tmp}METADATA/Reference_Sequences/{reference}/kraken2/{reftype}/opts.k2d",
        kraktaxo="{tmp}METADATA/Reference_Sequences/{reference}/kraken2/{reftype}/taxo.k2d",
        brakdb="{tmp}METADATA/Reference_Sequences/{reference}/kraken2/{reftype}/database1451mers.kraken",
        brakdist="{tmp}METADATA/Reference_Sequences/{reference}/kraken2/{reftype}/database1451mers.kmer_distrib"
    log:
        "{tmp}METADATA/Reference_Sequences/{reference}/kraken2/MeBaPiNa_{reftype}.log"
    benchmark:
        "{tmp}METADATA/Reference_Sequences/{reference}/kraken2/MeBaPiNa_{reftype}.benchmark.tsv"
    conda:
        "../envs/kraken2.yml"
    threads:
        4
    params:
        "35" ## k-mer length
    shell:
        "out_dir={output.krakdb}; out_dir=\"${{out_dir/database.kraken/}}\" > {log} 2>&1; "
        "kraken2-build --threads {threads} " ## unfortunately --kmer-len {params} is ignored by when --special is used
        "--special {wildcards.reference} --db ${{out_dir}} " ## reference can be one of "greengenes", "silva", "rdp"
        ">> {log} 2>&1; "
        "bracken-build -t {threads} -k {params} -l 1451 -d ${{out_dir}} " ## 1451 ismedian read length after filtering in 20191007_1559_MN31344_FAK76605_2bf006ff
        ">> {log} 2>&1; "
        "rm -rf ${{out_dir}}/data "
        "kraken2-build --clean --db ${{out_dir}} >> {log} 2>&1"

ruleorder: building_database_fromreffiles > building_database_alone

## MINIMA2 INDEX ##
###################

rule indexing_reference:
    input:
        "{tmp}METADATA/Reference_Sequences/{reference}/reference.fasta"
    output:
        "{tmp}METADATA/Reference_Sequences/{reference}/reference.mmi"
    log:
        "{tmp}METADATA/Reference_Sequences/{reference}/MeBaPiNa_indexing.log"
    benchmark:
        "{tmp}METADATA/Reference_Sequences/{reference}/MeBaPiNa_indexing.benchmark.tsv"
    params:
        "-x map-ont" ## naopore specific
    conda:
        "../envs/minimap2.yml"
    threads:
        2
    shell:
        "minimap2 -t {threads} {params} -d {output} {input} > {log} 2>&1"
