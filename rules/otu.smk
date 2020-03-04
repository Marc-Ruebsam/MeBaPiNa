###############
## IMPORTING ##
###############

rule generate_manifestfile:
    input:
        "{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/filtered.fastq"
    output:
        "{tmp}01_processed_data/03_otu_picking/{run}/{barc}/samplemanifest.tsv"
    run:
        ## import
        import re
        ## open file for writing
        wo = open(output[0], "w")
        
        ## write header
        null = wo.write("sample-id\tabsolute-filepath\n")
        ## write samples
        # for barc in SAMPLES:
        #     ## get input_file
        #     coorect_input_file = list(filter(lambda input_file: barc in input_file, input))[0]
        #     wo.write(SAMPLES[barc] + "\t" + "$PWD/" + coorect_input_file + "\n")
        wo.write(SAMPLES[wildcards.barc] + "\t" + "$PWD/" + input[0] + "\n")
        
        ## close file
        wo.close()

rule q2import_filtered:
    input:
        "{tmp}01_processed_data/03_otu_picking/{run}/{barc}/samplemanifest.tsv"
    output:
        "{tmp}01_processed_data/03_otu_picking/{run}/{barc}/filtered.qza"
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
        derepseq="{tmp}01_processed_data/03_otu_picking/{run}/{barc}/filtered_derep.qza",
        dereptable="{tmp}01_processed_data/03_otu_picking/{run}/{barc}/filtered_derep_table.qza"
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
        q2ref="{tmp}.tmp/Reference_Sequences/{reference}/qiime/reference.qza"
    output:
        otuseq="{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/cluster_seq.qza",
        otutable="{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/cluster_table.qza",
        otunewref="{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/cluster_newref.qza" ## The new reference sequences. This can be used for subsequent runs of open-reference clustering for consistent definitions of features across open-reference feature tables.
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

# qiime tools export --input-path cluster_table.qza --output-path .
# biom summarize-table -i feature-table.biom

# qiime tools import \
#   --type 'FeatureData[Sequence]' \
#   --input-path 85_otus.fasta \
#   --output-path 85_otus.qza
# 
# qiime tools import \
#   --type 'FeatureData[Taxonomy]' \
#   --input-format HeaderlessTSVTaxonomyFormat \
#   --input-path 85_otu_taxonomy.txt \
#   --output-path ref-taxonomy.qza

# qiime feature-classifier extract-reads \
#   --i-sequences 85_otus.qza \
#   --p-f-primer GTGCCAGCMGCCGCGGTAA \
#   --p-r-primer GGACTACHVGGGTWTCTAAT \
#   --p-trunc-len 120 \
#   --p-min-length 100 \
#   --p-max-length 400 \
#   --o-reads ref-seqs.qza

# qiime feature-classifier fit-classifier-naive-bayes \
#   --i-reference-reads ref-seqs.qza \
#   --i-reference-taxonomy ref-taxonomy.qza \
#   --o-classifier classifier.qza

# qiime feature-classifier classify-sklearn \
#   --i-classifier classifier.qza \
#   --i-reads rep-seqs.qza \
#   --o-classification taxonomy.qza
# 
# qiime metadata tabulate \
#   --m-input-file taxonomy.qza \
#   --o-visualization taxonomy.qzv
