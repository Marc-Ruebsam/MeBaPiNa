## KRAKEN2 ##

rule building_database:
    output:
        directory("METADATA/Reference_Sequences/kraken2_{reference}_132")
    log:
        "METADATA/Reference_Sequences/kraken2_{reference}_132/MeBaPiNa_building.log"
    benchmark:
        "METADATA/Reference_Sequences/kraken2_{reference}_132/MeBaPiNa_building.benchmark.tsv"
    conda:
        "../envs/kraken2.yml"
    threads:
        8
    params:
        "-k 35", ## k-mer length
        "-l 1451" ## median read length after filtering in 20191007_1559_MN31344_FAK76605_2bf006ff
    shell:
        # "kraken2-build --threads {threads} --download-taxonomy --skip-maps --db {output}; "
        # "kraken2-build --threads {threads} --download-library bacteria --no-masking --db {output}; "
        # "kraken2-build --threads {threads} --build --db {output}; "
        "kraken2-build --threads {threads} --special {wildcards.reference} --db {output} > {log} 2>&1; " ## reference can be one of "greengenes", "silva", "rdp"
        "bracken-build -d {output} -t {threads} {params} >> {log} 2>&1"
        # "kraken2-build --clean --db {output} >> {log} 2>&1"

rule kmer_mapping_filtered:
    input:
        fastq="01_processed_data/02_trimming_filtering/{run}/{barc}/filtered.fastq", 
        target=expand("METADATA/Reference_Sequences/kraken2_{reference}_132", reference=config["align"]["reference"])
    output:
        report="01_processed_data/03_kmer_mapping/{run}/{barc}/filtered.kreport2",
        output="01_processed_data/03_kmer_mapping/{run}/{barc}/filtered.kraken2"
    log:
        "01_processed_data/03_kmer_mapping/{run}/{barc}/MeBaPiNa_kmer_mapped.log"
    benchmark:
        "01_processed_data/03_kmer_mapping/{run}/{barc}/MeBaPiNa_kmer_mapped.benchmark.tsv"
    conda:
        "../envs/kraken2.yml"
    threads:
        8
    params:
        "--confidence 0.0 " ## how many of the k-mers have to map to a reference to be assigned (higher taxonomies accumulate the counts of lower ones)
    shell:
        "kraken2 --threads {threads} {params}"
        "--db {input.target} "
        # "--use-names " ## output contains scientific names and taxids
        # "--unclassified-out {output.unclas} "
        # "--classified-out {output.clas} "
        "--output {output.output} " ## information per sequence
        "--report {output.report} " ## information per taxon
        "{input.fastq} > {log} 2>&1"

rule ranking_taxonomy:
    input:
        report="01_processed_data/03_kmer_mapping/{run}/{barc}/filtered.kreport2",
        output="01_processed_data/03_kmer_mapping/{run}/{barc}/filtered.kraken2",
        target=expand("METADATA/Reference_Sequences/kraken2_{reference}_132", reference=config["align"]["reference"])
    output:
        rankS="02_analysis_results/03_kmer_mapping/{run}/bracken/{barc}/Species.bracken",
        tableS="02_analysis_results/03_kmer_mapping/{run}/bracken/{barc}/Species.breport",
        rankG="02_analysis_results/03_kmer_mapping/{run}/bracken/{barc}/Genus.bracken",
        tableG="02_analysis_results/03_kmer_mapping/{run}/bracken/{barc}/Genus.breport",
        rankF="02_analysis_results/03_kmer_mapping/{run}/bracken/{barc}/Family.bracken",
        tableF="02_analysis_results/03_kmer_mapping/{run}/bracken/{barc}/Family.breport"
    log:
        "02_analysis_results/03_kmer_mapping/{run}/bracken/{barc}/MeBaPiNa_ranking.log"
    benchmark:
        "02_analysis_results/03_kmer_mapping/{run}/bracken/{barc}/MeBaPiNa_ranking.benchmark.tsv"
    conda:
        "../envs/kraken2.yml"
    params:
        "-t 10", ## [Default = 10]:: specifies the minimum number of reads required for a classification at the specified rank. Any classifications with less than the specified threshold will not receive additional reads from higher taxonomy levels when distributing reads for abundance estimation.
        "-r 1451" ## median read length after filtering in 20191007_1559_MN31344_FAK76605_2bf006ff
    shell:
        "bracken {params} "
        "-l S " ## [Default = 'S', Options = 'D','P','C','O','F','G','S']:: specifies the taxonomic rank to analyze. Each classification at this specified rank will receive an estimated number of reads belonging to that rank after abundance estimation.
        "-d {input.target} "
        "-i {input.report} "
        "-o {output.rankS} > {log} 2>&1; "
        "mv 01_processed_data/{wildcards.run}/kmer_mapped/{wildcards.barc}/{wildcards.barc}_bracken.kreport2 {output.tableS}; "
        "bracken {params} "
        "-l G "
        "-d {input.target} "
        "-i {input.report} "
        "-o {output.rankG} >> {log} 2>&1; "
        "mv 01_processed_data/{wildcards.run}/kmer_mapped/{wildcards.barc}/{wildcards.barc}_bracken.kreport2 {output.tableG}; "
        "bracken {params} "
        "-l F "
        "-d {input.target} "
        "-i {input.report} "
        "-o {output.rankF} >> {log} 2>&1; "
        "mv 01_processed_data/{wildcards.run}/kmer_mapped/{wildcards.barc}/{wildcards.barc}_bracken.kreport2 {output.tableF} "
