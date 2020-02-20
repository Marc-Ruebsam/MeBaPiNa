## KRAKEN2 ##

rule kmer_mapping_filtered:
    input:
        fastq="{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/filtered.fastq", 
        target="{tmp}METADATA/Reference_Sequences/{reference}/kraken2/{reftype}"
    output:
        report="{tmp}01_processed_data/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/filtered.kreport2",
        output="{tmp}01_processed_data/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/filtered.kraken2"
    log:
        "{tmp}01_processed_data/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_kmer_mapped.log"
    benchmark:
        "{tmp}01_processed_data/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_kmer_mapped.benchmark.tsv"
    conda:
        "../envs/kraken2.yml"
    threads:
        8
    params:
        "--confidence 0.0" ## how many of the k-mers have to map to a reference to be assigned (higher taxonomies accumulate the counts of lower ones)
    shell:
        "kraken2 --threads {threads} {params} "
        "--db {input.target} "
        "--output {output.output} " ## information per sequence
        "--report {output.report} " ## information per taxon
        "{input.fastq} > {log} 2>&1"

rule ranking_taxonomy:
    input:
        report="{tmp}01_processed_data/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/filtered.kreport2",
        output="{tmp}01_processed_data/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/filtered.kraken2",
        target="{tmp}METADATA/Reference_Sequences/{reference}/kraken2/{reftype}"
    output:
        rankS="{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/bracken/{reference}_{reftype}/Species.bracken",
        tableS="{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/bracken/{reference}_{reftype}/Species.kreport2",
        rankG="{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/bracken/{reference}_{reftype}/Genus.bracken",
        tableG="{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/bracken/{reference}_{reftype}/Genus.kreport2",
        rankF="{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/bracken/{reference}_{reftype}/Family.bracken",
        tableF="{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/bracken/{reference}_{reftype}/Family.kreport2"
    log:
        "{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/bracken/{reference}_{reftype}/MeBaPiNa_ranking.log"
    benchmark:
        "{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/bracken/{reference}_{reftype}/MeBaPiNa_ranking.benchmark.tsv"
    conda:
        "../envs/kraken2.yml"
    params:
        "-t 3", ## [Default = 10]:: specifies the minimum number of reads required for a classification at the specified rank. Any classifications with less than the specified threshold will not receive additional reads from higher taxonomy levels when distributing reads for abundance estimation.
        "-r 1451" ## median read length after filtering in 20191007_1559_MN31344_FAK76605_2bf006ff
    shell:
        "bracken {params} "
        "-l S " ## [Default = 'S', Options = 'D','P','C','O','F','G','S']:: specifies the taxonomic rank to analyze. Each classification at this specified rank will receive an estimated number of reads belonging to that rank after abundance estimation.
        "-d {input.target} "
        "-i {input.report} "
        "-o {output.rankS} > {log} 2>&1; "
        "mv {wildcards.tmp}01_processed_data/03_kmer_mapping/{wildcards.run}/{wildcards.barc}/{wildcards.reference}_{wildcards.reftype}/filtered_bracken.kreport2 {output.tableS}; "
        "bracken {params} "
        "-l G "
        "-d {input.target} "
        "-i {input.report} "
        "-o {output.rankG} >> {log} 2>&1; "
        "mv {wildcards.tmp}01_processed_data/03_kmer_mapping/{wildcards.run}/{wildcards.barc}/{wildcards.reference}_{wildcards.reftype}/filtered_bracken.kreport2 {output.tableG}; "
        "bracken {params} "
        "-l F "
        "-d {input.target} "
        "-i {input.report} "
        "-o {output.rankF} >> {log} 2>&1; "
        "mv {wildcards.tmp}01_processed_data/03_kmer_mapping/{wildcards.run}/{wildcards.barc}/{wildcards.reference}_{wildcards.reftype}/filtered_bracken.kreport2 {output.tableF} "
