###############
## IMPORTING ##
###############

rule generate_manifestfile:
    input:
        "{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/filtered.fastq"
    output:
        temp("{tmp}01_processed_data/03_otu_picking/{run}/{barc}/samplemanifest.tsv")
    run:
        ## import
        import re
        ## open file for writing
        wo = open(output[0], "w")
        ## write header
        null = wo.write("sample-id\tabsolute-filepath\n")
        ## write samples
        wo.write(SAMPLES[wildcards.barc] + "\t" + "$PWD/" + input[0] + "\n")
        ## close file
        wo.close()

rule q2import_filtered:
    input:
        "{tmp}01_processed_data/03_otu_picking/{run}/{barc}/samplemanifest.tsv"
    output:
        temp("{tmp}01_processed_data/03_otu_picking/{run}/{barc}/filtered.qza")
    log:
        "{tmp}01_processed_data/03_otu_picking/{run}/{barc}/MeBaPiNa_q2import_filtered.log"
    benchmark:
        "{tmp}01_processed_data/03_otu_picking/{run}/{barc}/MeBaPiNa_q2import_filtered.benchmark.tsv"
    conda:
        "../envs/qiime2.yml"
    threads:
        1
    params:
        "--type 'SampleData[SequencesWithQuality]'",
        "--input-format SingleEndFastqManifestPhred33V2"
    shell:
        "qiime tools import {params} "
        "--input-path {input} "
        "--output-path {output} > {log} 2>&1"

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
    conda:
        "../envs/qiime2.yml"
    threads:
        1
    params:
        "--verbose"
    shell:
        "qiime vsearch dereplicate-sequences {params} "
        "--i-sequences {input} "
        "--o-dereplicated-table {output.dereptable} "
        "--o-dereplicated-sequences {output.derepseq} "
        "> {log} 2>&1"

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
    conda:
        "../envs/qiime2.yml"
    params:
        "--p-perc-identity " + config['filtering']['min_readidentity'], ## The percent identity at which clustering should be performed.
        "--verbose"
    threads:
        16 ## job needs some tmp storage space, to many parallel jobs might break stuff!?
    shell:
        "qiime vsearch cluster-features-open-reference --p-threads {threads} {params} "
        "--i-sequences {input.derepseq} --i-table {input.dereptable} "
        "--i-reference-sequences {input.q2ref} " ## FeatureData[Sequence]
        "--o-clustered-sequences {output.otuseq} "
        "--o-clustered-table {output.otutable} "
        "--o-new-reference-sequences {output.otunewref} "
        "> {log} 2>&1"

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
    conda:
        "../envs/qiime2.yml"
    params:
        "--verbose"
    threads:
        1
    shell:
        "qiime vsearch uchime-denovo {params} "
        "--i-table {input.otutable} "
        "--i-sequences {input.otuseq} "
        "--o-nonchimeras {output.nonchim} "
        "--o-chimeras {output.chimera} "
        "--o-stats {output.stats} "
        "> {log} 2>&1"

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
    conda:
        "../envs/qiime2.yml"
    params:
        "--p-min-frequency " + config['filtering']['min_featurereads'] ## The minimum total frequency that a feature must have to be retained.
    threads:
        1
    shell:
        "qiime feature-table filter-features " ## exclude all chimeras from feature count table
        "--i-table {input.otutable} "
        "--m-metadata-file {input.chimera} "
        "--p-exclude-ids "
        "--o-filtered-table {output.nochimtable} "
        "--verbose {params} " ## params excludes clusters with low counts
        "> {log} 2>&1; "
        "qiime feature-table filter-seqs " ## retain only sequences from ids retained above
        "--i-data {input.otuseq} "
        "--i-table {output.nochimtable} "
        "--p-no-exclude-ids "
        "--o-filtered-data {output.nochimseq} "
        "--verbose "
        ">> {log} 2>&1"

################
## CONVERSION ##
################

rule convert_q2ftable:
    input:
        "{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/{ftable}_ftable.qza"
    output:
        ftable=temp("{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/{ftable}_ftable/feature-table.tsv"),
        ftablebiom=temp("{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/{ftable}_ftable/temp.biom") #!# separated coversion rules to prevent rerunning of rules rereplicate_q2filte, ... when rerunning rule stat_otu_feature
    log:
        "{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/{ftable}_ftable/MeBaPiNa_convert_q2ftable.log"
    benchmark:
        "{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/{ftable}_ftable/MeBaPiNa_convert_q2ftable.benchmark.tsv"
    conda:
        "../envs/qiime2.yml"
    shell:
        "qiime tools export --input-path {input} --output-path {output.ftablebiom} > {log} 2>&1; "
        "biom convert --input-fp {output.ftablebiom} --output-fp {output.ftable} --to-tsv >> {log} 2>&1"

rule convert_q2ftable_biom:
    input:
        "{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/{ftable}_ftable.qza"
    output:
        temp("{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/{ftable}_ftable/feature-table.biom") #!# separated coversion rules to prevent rerunning of rules rereplicate_q2filte, ... when rerunning rule stat_otu_feature
    log:
        "{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/{ftable}_ftable/MeBaPiNa_convert_q2ftable_biom.log"
    benchmark:
        "{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/{ftable}_ftable/MeBaPiNa_convert_q2ftable_biom.benchmark.tsv"
    conda:
        "../envs/qiime2.yml"
    shell:
        "qiime tools export --input-path {input} --output-path {output} > {log} 2>&1"

rule convert_q2centseq:
    input:
        "{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/{centseq}_centseq.qza"
    output:
        centseq=temp("{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/{centseq}_centseq/dna-sequences.fasta")
    log:
        "{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/{ftable}_ftable/MeBaPiNa_convert_q2centseq.log"
    benchmark:
        "{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/{ftable}_ftable/MeBaPiNa_convert_q2centseq.benchmark.tsv"
    conda:
        "../envs/qiime2.yml"
    shell:
        "qiime tools export --input-path {input} --output-path {output} > {log} 2>&1"

##########################
## TAXONOMIC ASSIGNMENT ##
##########################

# rule q2filter_classify:
#     input:
#         nochimseq="{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/filt_centseq.qza",
#         classifier="{tmp}METADATA/Reference_Sequences/silva/qiime/{reftype}/classifyer.qza" ## not available
#     output:
#         "{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}_{reftype}/counttax.qza"
#     log:
#         "{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_q2filter_classify.log"
#     benchmark:
#         "{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_q2filter_classify.benchmark.tsv"
#     conda:
#         "../envs/qiime2.yml"
#     params:
#         "--p-confidence " + config["filtering"]["min_confidence"] ## Confidence threshold for limiting taxonomic depth. Set to "disable" to disable confidence calculation, or 0 to calculate confidence but not apply it to limit the taxonomic depth of the assignments. [default: 0.7]
#     threads:
#         8
#     shell:
#         "qiime feature-classifier classify-sklearn "
#         "--i-reads {input.nochimseq} "
#         "--i-classifier {input.classifier} "
#         "--o-classification {output} "
#         "--p-n-jobs {threads} "
#         "--verbose {params} >> {log} 2>&1"

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
    threads:
        1
    shell:
        "awk 'FNR==NR&&!/^#/{{featcount[$1]=$2}}; " ## get number of reads per featureID from ftable
        "FNR!=NR&&/^>/{{cur_id=substr($1,2)}}; " ## get current featureIF in fasta
        "FNR!=NR&&!/^>/{{for(i=0;i<featcount[cur_id];i++){{print \">\"cur_id\" \"i\"\\n\"$0}}}}' " ## for current featureID: print featureID and sequene times the number of reads
        "{input.ftable} {input.centseq} > {output.rerep} 2> {log}"

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
    conda:
        "../envs/kraken2.yml"
    threads:
        8
    params:
        "--confidence " + config["filtering"]["min_confidence"] ## how many of the k-mers have to map to a reference to be assigned (higher taxonomies accumulate the counts of lower ones)
    shell:
        "target={input.krakdb}; target=\"${{target/database.kraken/}}\" >> {log} 2>&1; "
        "kraken2 --threads {threads} {params} "
        "--db ${{target}} "
        "--output {output.output} " ## information per sequence
        "--report {output.report} " ## information per taxon
        "{input.rerep} >> {log} 2>&1"

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
    conda:
        "../envs/python.yml"
    script:
        "../scripts/convert_kreport.py"
