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
        "{tmp}01_processed_data/01_basecalling/{run}/sequencing_summary/MeBaPiNa_sorting_seqsum_barc.log"
    benchmark:
        "{tmp}01_processed_data/01_basecalling/{run}/sequencing_summary/MeBaPiNa_sorting_seqsum_barc.benchmark.tsv"
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
        "{tmp}01_processed_data/01_basecalling/{run}/sequencing_summary/MeBaPiNa_splitting_seqsum_barc.log"
    benchmark:
        "{tmp}01_processed_data/01_basecalling/{run}/sequencing_summary/MeBaPiNa_splitting_seqsum_barc.benchmark.tsv"
    shell:
        "mkdir {output[0]}"

## DOWNSAMPLING ##
##################

rule downsampling_seqsum: #!#
    input:
        "{tmp}01_processed_data/01_basecalling/{run}/sequencing_summary/sequencing_summary_sorted.txt"
    output:
        temp("{tmp}01_processed_data/01_basecalling/{run}/sequencing_summary/sequencing_summary_downsampled.txt")
    log:
        "{tmp}01_processed_data/01_basecalling/{run}/sequencing_summary/MeBaPiNa_downsampling_seqsum.log"
    benchmark:
        "{tmp}01_processed_data/01_basecalling/{run}/sequencing_summary/MeBaPiNa_downsampling_seqsum.benchmark.tsv"
    shell:
        "touch {output[0]}"

##############
## SAMTOOLS ##
##############

rule samtools_faidx:
    input:
        "{fasta}"
    output:
        temp("{fasta}.fai")
    conda:
        "../envs/samtools.yml"
    shell:
        "samtools faidx {input}"

rule samtools_depth:
    input:
        bam="{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/filteredsorted.bam",
        reference="{tmp}METADATA/Reference_Sequences/{reference}/reference.fasta"
    output:
        temp("{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/refseq_coverage.tsv")
    conda:
        "../envs/samtools.yml"
    params:
        "-a",
        "-d 0", ## maximum coverage depth [8000]. If 0, depth is set to the maximum
        "-l " + config['filtering']['len_min'] #!# should be redundant
    shell:
        "samtools depth {params} "
        "--reference {input.reference} {input.bam} "
        "> {output}"

################
## REFERENCES ##
################

## SILVA ##
###########

rule download_reffiles:
    output:
        fasta_align=temp("{tmp}METADATA/Reference_Sequences/silva/reference_aligned.fasta"), ## uses the non-redundant (99% clustering) file
        slvmap=temp("{tmp}METADATA/Reference_Sequences/silva/slvmap.txt"),
        taxlist=temp("{tmp}METADATA/Reference_Sequences/silva/taxlist.txt")
    log:
        "{tmp}METADATA/Reference_Sequences/silva/MeBaPiNa_download_reffiles.log"
    benchmark:
        "{tmp}METADATA/Reference_Sequences/silva/MeBaPiNa_download_reffiles.benchmark.tsv"
    shell:
        "bash {wildcards.tmp}Pipeline/MeBaPiNa/scripts/download_silva.sh {wildcards.tmp}METADATA/Reference_Sequences/silva > {log} 2>&1"

rule construct_refseq:
    input:
        refseq="{tmp}METADATA/Reference_Sequences/silva/reference_aligned.fasta",
        primers="{tmp}METADATA/Reference_Sequences/primers/" + config['primers'] + "/primers.fasta"
    output:
        refseq="{tmp}METADATA/Reference_Sequences/silva/reference.fasta",
        dups=temp("{tmp}METADATA/Reference_Sequences/silva/reference.dups")
    log:
        "{tmp}METADATA/Reference_Sequences/silva/MeBaPiNa_construct_refseq.log"
    benchmark:
        "{tmp}METADATA/Reference_Sequences/silva/MeBaPiNa_construct_refseq.benchmark.tsv"
    conda:
        "../envs/qiime2.yml"
    threads:
        2
    params:
        "--fastq_minlen " + config["filtering"]["len_min"], ## minimum length of amplicon
        "--fastq_maxlen " + config["filtering"]["len_max"]
    shell:
        ## convert sequences to DNA with gaps
        "python {wildcards.tmp}Pipeline/MeBaPiNa/scripts/make_SILVA_db/convert_rna_to_dna.py "
        "--convert_to_gap "
        "--input_fasta {input.refseq} "
        "--output_fasta {output.refseq}_dna.fasta "
        "> {log} 2>&1; "
        ## extract first 500 sequences
        "head -n 500 {output.refseq}_dna.fasta > {output.refseq}_head.fasta 2>> {log}; "
        ## align primers to first 500 sequences (aka add primer sequences to multiple alignment)
        "mafft --thread {threads} "
        "--addfragments {input.primers} "
        "--mapout {output.refseq}_head.fasta  > /dev/null 2>> {log}; "
        "rm {output.refseq}_head.fasta; "
        ## trim ALL sequences to region of primer match
        "python {wildcards.tmp}Pipeline/MeBaPiNa/scripts/make_SILVA_db/extract_alignment_region.py "
        "--input_alignment {output.refseq}_dna.fasta "
        "--output_alignment {output.refseq}_trim.fasta "
        "--start_position $(awk 'BEGIN{{prnt=-1}}; />/{{prnt++}}; prnt==0{{print $3+1}}' {input.primers}.map | tail -n 1) "
        "--end_position $(awk 'BEGIN{{prnt=-3}}; />/{{prnt++}}; prnt==0{{print $3-1}}' {input.primers}.map | head -n 3 | tail -n 1) "
        ">> {log} 2>&1; "
        "rm {output.refseq}_dna.fasta; "
        "rm {input.primers}.map; "
        ## remove gaps
        "python {wildcards.tmp}Pipeline/MeBaPiNa/scripts/make_SILVA_db/degap_fasta.py "
        "--input_fasta {output.refseq}_trim.fasta "
        "--output_fasta {output.refseq}_degap.fasta "
        ">> {log} 2>&1; "
        "rm {output.refseq}_trim.fasta; "
        ## exclude sequences above or below the thresholds
        "vsearch --fastx_filter "
        "{output.refseq}_degap.fasta "
        "--fastaout {output.refseq}_filt.fasta {params} "
        ">> {log} 2>&1; "
        "rm {output.refseq}_degap.fasta; "
        ## extract ids and group duplicates (first id is the one kept after deduplication)
        "awk 'BEGIN{{init=0}}; "
        "/^>/&&init==0{{ accID = $0; init++; next }}; "
        "/^>/{{ accID = substr(accID,2); print accID\"\\t\"sqnc; "
        "accID = $0; sqnc = \"\"; next }}; "
        "{{ sqnc = sqnc$0 }}' {output.refseq}_filt.fasta | "
        "awk '{{ ids[$2] = ids[$2]$1\";\" }}; "
        "END{{ for( sqnc in ids ){{ "
        "print ids[sqnc] }} }}' > {output.dups} "
        "2>> {log}; "
        ## remove replicates (deduplication)
        "vsearch --derep_fulllength "
        "{output.refseq}_filt.fasta "
        "--output {output.refseq} "
        ">> {log} 2>&1; "
        "rm {output.refseq}_filt.fasta"

rule construct_reftax:
    input:
        taxlist="{tmp}METADATA/Reference_Sequences/silva/taxlist.txt",
        slvmap="{tmp}METADATA/Reference_Sequences/silva/slvmap.txt",
        dups="{tmp}METADATA/Reference_Sequences/silva/reference.dups"
    output:
        kraknames_S="{tmp}METADATA/Reference_Sequences/silva/kraken2/Species/taxonomy/names.dmp", ## removed after reference build
        kraknodes_S="{tmp}METADATA/Reference_Sequences/silva/kraken2/Species/taxonomy/nodes.dmp", ## removed after reference build
        krakseq2tax_S="{tmp}METADATA/Reference_Sequences/silva/kraken2/Species/seqid2taxid.map", ## removed after reference build
        kraknames_G="{tmp}METADATA/Reference_Sequences/silva/kraken2/Genus/taxonomy/names.dmp", ## removed after reference build
        kraknodes_G="{tmp}METADATA/Reference_Sequences/silva/kraken2/Genus/taxonomy/nodes.dmp", ## removed after reference build
        krakseq2tax_G="{tmp}METADATA/Reference_Sequences/silva/kraken2/Genus/seqid2taxid.map", ## removed after reference build
        kronataxtab_S="{tmp}METADATA/Reference_Sequences/silva/krona/Species/taxonomy.tab",
        kronataxlist_S="{tmp}METADATA/Reference_Sequences/silva/krona/Species/taxlist.txt",
        kronaseq2tax_S="{tmp}METADATA/Reference_Sequences/silva/krona/Species/seqid2taxid.map",
        kronataxtab_G="{tmp}METADATA/Reference_Sequences/silva/krona/Genus/taxonomy.tab",
        kronataxlist_G="{tmp}METADATA/Reference_Sequences/silva/krona/Genus/taxlist.txt",
        kronaseq2tax_G="{tmp}METADATA/Reference_Sequences/silva/krona/Genus/seqid2taxid.map",
        qiimetax_S=temp("{tmp}METADATA/Reference_Sequences/silva/qiime/Species/taxonomy.tsv"),
        qiimetax_G=temp("{tmp}METADATA/Reference_Sequences/silva/qiime/Genus/taxonomy.tsv")
    log:
        "{tmp}METADATA/Reference_Sequences/silva/MeBaPiNa_construct_reftax.log"
    benchmark:
        "{tmp}METADATA/Reference_Sequences/silva/MeBaPiNa_construct_reftax.benchmark.tsv"
    conda:
        "../envs/python.yml"
    script:
        "../scripts/construct_reffiles.py"

## QIIME REFERENCE ##
#####################

rule q2import_reftax:
    input:
        "{tmp}METADATA/Reference_Sequences/{reference}/qiime/{reftype}/taxonomy.tsv"
    output:
        "{tmp}METADATA/Reference_Sequences/{reference}/qiime/{reftype}/taxonomy.qza"
    log:
        "{tmp}METADATA/Reference_Sequences/{reference}/qiime/{reftype}/MeBaPiNa_q2import_reftax.log"
    benchmark:
        "{tmp}METADATA/Reference_Sequences/{reference}/qiime/{reftype}/MeBaPiNa_q2import_reftax.benchmark.tsv"
    conda:
        "../envs/qiime2.yml"
    threads:
        1
    shell:
        "qiime tools import {params} "
        "--type FeatureData[Taxonomy] "
        "--input-format HeaderlessTSVTaxonomyFormat "
        "--input-path {input} "
        "--output-path {output} "
        "> {log} 2>&1"

rule q2import_refseq:
    input:
        "{tmp}METADATA/Reference_Sequences/{reference}/reference.fasta"
    output:
        "{tmp}METADATA/Reference_Sequences/{reference}/qiime/reference.qza"
    log:
        "{tmp}METADATA/Reference_Sequences/{reference}/qiime/MeBaPiNa_q2import_refseq.log"
    benchmark:
        "{tmp}METADATA/Reference_Sequences/{reference}/qiime/MeBaPiNa_q2import_refseq.benchmark.tsv"
    conda:
        "../envs/qiime2.yml"
    threads:
        1
    shell:
        "qiime tools import {params} "
        "--type FeatureData[Sequence] "
        "--input-path {input} "
        "--output-path {output} "
        "> {log} 2>&1"

rule q2train_classifyer: #!# failed repeatedly because of high memory requirements
    input:
        reftax="{tmp}METADATA/Reference_Sequences/{reference}/qiime/{reftype}/taxonomy.qza",
        refseq="{tmp}METADATA/Reference_Sequences/{reference}/qiime/reference.qza"
    output:
        "{tmp}METADATA/Reference_Sequences/{reference}/qiime/{reftype}/classifyer.qza"
    log:
        "{tmp}METADATA/Reference_Sequences/{reference}/qiime/{reftype}/MeBaPiNa_q2train_classifyer.log"
    benchmark:
        "{tmp}METADATA/Reference_Sequences/{reference}/qiime/{reftype}/MeBaPiNa_q2train_classifyer.benchmark.tsv"
    conda:
        "../envs/qiime2.yml"
    threads:
        1
    params:
        "--verbose"
    shell:
        "qiime feature-classifier fit-classifier-naive-bayes {params} "
        "--i-reference-reads {input.refseq} "
        "--i-reference-taxonomy {input.reftax} "
        "--o-classifier {output} "
        "> {log} 2>&1"

## MINIMAP2 INDEX ##
###################

rule indexing_reference:
    input:
        "{tmp}METADATA/Reference_Sequences/{reference}/reference.fasta"
    output:
        "{tmp}METADATA/Reference_Sequences/{reference}/reference.mmi"
    log:
        "{tmp}METADATA/Reference_Sequences/{reference}/MeBaPiNa_indexing_reference.log"
    benchmark:
        "{tmp}METADATA/Reference_Sequences/{reference}/MeBaPiNa_indexing_reference.benchmark.tsv"
    params:
        "-x map-ont" ## naopore specific
    conda:
        "../envs/minimap2.yml"
    threads:
        2
    shell:
        "minimap2 -t {threads} {params} -d {output} {input} > {log} 2>&1"

## KRAKEN2 DATABASE ##
######################

rule building_database:
    input:
        fasta="{tmp}METADATA/Reference_Sequences/silva/reference.fasta",
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
        "{tmp}METADATA/Reference_Sequences/{reference}/kraken2/MeBaPiNa_building_{reftype}_database.log"
    benchmark:
        "{tmp}METADATA/Reference_Sequences/{reference}/kraken2/MeBaPiNa_building_{reftype}_database.benchmark.tsv"
    conda:
        "../envs/kraken2.yml"
    threads:
        6
    params:
        "35" ## k-mer length
    shell:
        "out_dir={output.krakdb}; out_dir=\"${{out_dir/database.kraken/}}\" > {log} 2>&1; "
        "mkdir -p \"${{out_dir}}library\" >> {log} 2>&1; "
        "cp \"{input.fasta}\" \"${{out_dir}}library/library.fna\" >> {log} 2>&1; "
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

# rule building_database_plain:
#     output:
#         krakdb="{tmp}METADATA/Reference_Sequences/{reference}/kraken2/{reftype}/database.kraken", ## reftype wildcard has no influence on this rule, but is used for silva Species database (rule building_database)
#         krakhash="{tmp}METADATA/Reference_Sequences/{reference}/kraken2/{reftype}/hash.k2d",
#         krakopts="{tmp}METADATA/Reference_Sequences/{reference}/kraken2/{reftype}/opts.k2d",
#         kraktaxo="{tmp}METADATA/Reference_Sequences/{reference}/kraken2/{reftype}/taxo.k2d",
#         brakdb="{tmp}METADATA/Reference_Sequences/{reference}/kraken2/{reftype}/database1451mers.kraken",
#         brakdist="{tmp}METADATA/Reference_Sequences/{reference}/kraken2/{reftype}/database1451mers.kmer_distrib"
#     log:
#         "{tmp}METADATA/Reference_Sequences/{reference}/kraken2/MeBaPiNa_building_{reftype}_database.log"
#     benchmark:
#         "{tmp}METADATA/Reference_Sequences/{reference}/kraken2/MeBaPiNa_building_{reftype}_database.benchmark.tsv"
#     conda:
#         "../envs/kraken2.yml"
#     threads:
#         4
#     params:
#         "35" ## k-mer length
#     shell:
#         "out_dir={output.krakdb}; out_dir=\"${{out_dir/database.kraken/}}\" > {log} 2>&1; "
#         "kraken2-build --threads {threads} " ## unfortunately --kmer-len {params} is ignored by when --special is used
#         "--special {wildcards.reference} --db ${{out_dir}} " ## reference can be one of "greengenes", "silva", "rdp"
#         ">> {log} 2>&1; "
#         "bracken-build -t {threads} -k {params} -l 1451 -d ${{out_dir}} " ## 1451 ismedian read length after filtering in 20191007_1559_MN31344_FAK76605_2bf006ff
#         ">> {log} 2>&1; "
#         "rm -rf ${{out_dir}}/data "
#         "kraken2-build --clean --db ${{out_dir}} >> {log} 2>&1"
# ruleorder: building_database > building_database_plain
