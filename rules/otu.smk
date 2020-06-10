###############
## IMPORTING ##
###############

rule generate_manifestfile:
    input:
        "{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/filtered.fastq"
    output:
        temp("{tmp}01_processed_data/03_otu_picking/{run}/{barc}/samplemanifest.tsv")
    shell:
        "touch {output[0]}"

rule q2import_filtered:
    input:
        "{tmp}01_processed_data/03_otu_picking/{run}/{barc}/samplemanifest.tsv"
    output:
        temp("{tmp}01_processed_data/03_otu_picking/{run}/{barc}/filtered.qza")
    log:
        "{tmp}01_processed_data/03_otu_picking/{run}/{barc}/MeBaPiNa_q2import_filtered.log"
    benchmark:
        "{tmp}01_processed_data/03_otu_picking/{run}/{barc}/MeBaPiNa_q2import_filtered.benchmark.tsv"
    shell:
        "touch {output[0]}; "
        "lg={log}; "
        "cat ${{lg/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} > {log}"

################
## CLUSTERING ##
################

rule q2derep_imported:
    input:
        "{tmp}01_processed_data/03_otu_picking/{run}/{barc}/filtered.qza"
    output:
        derepseq=temp("{tmp}01_processed_data/03_otu_picking/{run}/{barc}/filtered_derep.qza"),
        dereptable=temp("{tmp}01_processed_data/03_otu_picking/{run}/{barc}/filtered_derep_table.qza")
    log:
        "{tmp}01_processed_data/03_otu_picking/{run}/{barc}/MeBaPiNa_q2derep_imported.log"
    benchmark:
        "{tmp}01_processed_data/03_otu_picking/{run}/{barc}/MeBaPiNa_q2derep_imported.benchmark.tsv"
    shell:
        "touch {output[0]}; "
        "touch {output[1]}; "
        "lg={log}; "
        "cat ${{lg/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} > {log}"

rule q2otupick:
    input:
        derepseq="{tmp}01_processed_data/03_otu_picking/{run}/{barc}/filtered_derep.qza",
        dereptable="{tmp}01_processed_data/03_otu_picking/{run}/{barc}/filtered_derep_table.qza",
        q2ref="{tmp}METADATA/Reference_Sequences/{reference}/qiime/reference.qza"
    output:
        otuseq="{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/cluster_centseq.qza", ## FeatureData[Sequence]
        otutable="{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/cluster_ftable.qza", ## FeatureTable[Frequency]
        otunewref="{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/cluster_newrefseq.qza" ## The new reference sequences. This can be used for subsequent runs of open-reference clustering for consistent definitions of features across open-reference feature tables.
    log:
        "{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/MeBaPiNa_q2otupick.log"
    benchmark:
        "{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/MeBaPiNa_q2otupick.benchmark.tsv"
    shell:
        "out1={output[0]}; "
        "mv ${{out1/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} ${{out1}}; "
        "out2={output[1]}; "
        "mv ${{out2/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} ${{out2}}; "
        "out3={output[2]}; "
        "mv ${{out3/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} ${{out3}}; "
        "lg={log}; "
        "cat ${{lg/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} > {log}"

################
## PROCESSING ##
################

rule q2uchime_otus:
    input:
        otuseq="{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/cluster_centseq.qza",
        otutable="{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/cluster_ftable.qza"
    output:
        nonchim=temp("{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/nochim_centseq.qza"),
        chimera=temp("{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/chim_centseq.qza"),
        stats=temp("{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/chim_stats.qza")
    log:
        "{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/MeBaPiNa_q2uchime_otus.log"
    benchmark:
        "{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/MeBaPiNa_q2uchime_otus.benchmark.tsv"
    shell:
        "touch {output[0]}; "
        "touch {output[1]}; "
        "touch {output[2]}; "
        "lg={log}; "
        "cat ${{lg/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} > {log}"

rule q2filter_uchime:
    input:
        otuseq="{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/cluster_centseq.qza",
        otutable="{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/cluster_ftable.qza",
        chimera="{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/chim_centseq.qza",
        nonchim="{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/nochim_centseq.qza"
    output:
        nochimtable="{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/filt_ftable.qza",
        nochimseq="{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/filt_centseq.qza"
    log:
        "{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/MeBaPiNa_q2filter_uchime.log"
    benchmark:
        "{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/MeBaPiNa_q2filter_uchime.benchmark.tsv"
    shell:
        "out1={output[0]}; "
        "mv ${{out1/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} ${{out1}}; "
        "out2={output[1]}; "
        "mv ${{out2/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} ${{out2}}; "
        "lg={log}; "
        "cat ${{lg/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} > {log}"

##########################
## TAXONOMIC ASSIGNMENT ##
##########################

rule convert_q2ftable:
    input:
        "{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/{ftable}_ftable.qza"
    output:
        ftable=temp("{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/{ftable}_ftable/feature-table.tsv"),
        ftablebiom=temp("{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/{ftable}_ftable/feature-table.biom")
    shell:
        "touch {output[0]}; "
        "touch {output[1]}"

rule convert_q2centseq:
    input:
        "{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/{centseq}_centseq.qza"
    output:
        centseq=temp("{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/{centseq}_centseq/dna-sequences.fasta")
    shell:
        "touch {output[0]}"

rule rereplicate_q2filter:
    input:
        ftable="{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/filt_ftable/feature-table.tsv",
        centseq="{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/filt_centseq/dna-sequences.fasta"
    output:
        rerep=temp("{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/filt_centseq/dna-sequences-rerep.fasta")
    log:
        "{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/MeBaPiNa_rereplicate_q2filter.log"
    benchmark:
        "{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/MeBaPiNa_rereplicate_q2filter.benchmark.tsv"
    shell:
        "touch {output[0]}; "
        "lg={log}; "
        "cat ${{lg/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} > {log}"

rule kmermap_q2rereplicate:
    input:
        rerep="{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/filt_centseq/dna-sequences-rerep.fasta",
        krakdb="{tmp}METADATA/Reference_Sequences/{reference}/kraken2/{reftype}/database.kraken"
    output:
        report="{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}_{reftype}/filtered.kreport2",
        output="{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}_{reftype}/filtered.kraken2"
    log:
        "{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_kmermap_q2rereplicate.log"
    benchmark:
        "{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_kmermap_q2rereplicate.benchmark.tsv"
    shell:
        "out1={output[0]}; "
        "mv ${{out1/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} ${{out1}}; "
        "out2={output[1]}; "
        "mv ${{out2/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} ${{out2}}; "
        "lg={log}; "
        "cat ${{lg/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} > {log}"

rule counttax_q2kmermap:
    input:
        kreport="{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}_{reftype}/filtered.kreport2",
        kronataxlist="{tmp}METADATA/Reference_Sequences/{reference}/krona/{reftype}/taxlist.txt"
    output:
        counttaxlist="{tmp}02_analysis_results/03_otu_picking/{run}/{barc}/{reference}_{reftype}/kmer.counttaxlist"
    log:
        "{tmp}02_analysis_results/03_otu_picking/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_counttax_q2kmermap.log"
    benchmark:
        "{tmp}02_analysis_results/03_otu_picking/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_counttax_q2kmermap.benchmark.tsv"
    shell:
        "out1={output[0]}; "
        "mv ${{out1/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} ${{out1}}; "
        "lg={log}; "
        "cat ${{lg/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} > {log}"
