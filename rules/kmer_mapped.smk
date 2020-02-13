## KRAKEN2 ##

rule kmer_mapping_filtered:
    input:
        fastq="01_processed_data/02_trimming_filtering/{run}/{barc}/filtered.fastq", 
        target="METADATA/Reference_Sequences/{reference}/kraken2/{reftype}"
    output:
        report="01_processed_data/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/filtered.kreport2",
        output="01_processed_data/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/filtered.kraken2"
    log:
        "01_processed_data/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_kmer_mapped.log"
    benchmark:
        "01_processed_data/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_kmer_mapped.benchmark.tsv"
    conda:
        "../envs/kraken2.yml"
    threads:
        8
    params:
        "--confidence 0.0 " ## how many of the k-mers have to map to a reference to be assigned (higher taxonomies accumulate the counts of lower ones)
    shell:
        "kraken2 --threads {threads} {params}"
        "--db {input.target} "
        "--output {output.output} " ## information per sequence
        "--report {output.report} " ## information per taxon
        "{input.fastq} > {log} 2>&1"

rule ranking_taxonomy:
    input:
        report="01_processed_data/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/filtered.kreport2",
        output="01_processed_data/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/filtered.kraken2",
        target=expand("METADATA/Reference_Sequences/kraken2_{reference}", reference = config["reference"]["source"])
    output:
        rankS="02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/bracken/Species.bracken",
        tableS="02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/bracken/Species.kreport2",
        rankG="02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/bracken/Genus.bracken",
        tableG="02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/bracken/Genus.kreport2",
        rankF="02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/bracken/Family.bracken",
        tableF="02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/bracken/Family.kreport2"
    log:
        "02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/bracken/MeBaPiNa_ranking.log"
    benchmark:
        "02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/bracken/MeBaPiNa_ranking.benchmark.tsv"
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
        "mv 01_processed_data/03_kmer_mapping/{wildcards.run}/{wildcards.barc}/{reference}_{reftype}/filtered_bracken.kreport2 {output.tableS}; "
        "bracken {params} "
        "-l G "
        "-d {input.target} "
        "-i {input.report} "
        "-o {output.rankG} >> {log} 2>&1; "
        "mv 01_processed_data/03_kmer_mapping/{wildcards.run}/{wildcards.barc}/{reference}_{reftype}/filtered_bracken.kreport2 {output.tableG}; "
        "bracken {params} "
        "-l F "
        "-d {input.target} "
        "-i {input.report} "
        "-o {output.rankF} >> {log} 2>&1; "
        "mv 01_processed_data/03_kmer_mapping/{wildcards.run}/{wildcards.barc}/{reference}_{reftype}/filtered_bracken.kreport2 {output.tableF} "
