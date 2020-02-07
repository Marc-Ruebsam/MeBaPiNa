## KRAKEN2 ##

rule building_database:
    output:
        directory("METADATA/Reference_Sequences/kraken2_{reference}")
    log:
        "METADATA/Reference_Sequences/kraken2_{reference}_132/MeBaPiNa_building.log"
    benchmark:
        "METADATA/Reference_Sequences/kraken2_{reference}_132/MeBaPiNa_building.benchmark.tsv"
    conda:
        "../envs/kraken2.yml"
    threads:
        8
    params:
        "35" ## k-mer length
    shell:
        # "kraken2-build --threads {threads} --download-taxonomy --skip-maps --db {output} > {log} 2>&1; "
        # "kraken2-build --threads {threads} --download-library bacteria --no-masking --db {output} >> {log} 2>&1; "
        # # "kraken2-build --threads {threads} --download-library archaea --no-masking --db {output} >> {log} 2>&1; "
        # "kraken2-build --threads {threads} --kmer-len {params} --build --db {output} >> {log} 2>&1; "
        # "kraken2-build --clean --db {output} >> {log} 2>&1"
        "kraken2-build --threads {threads} " ## unfortunately --kmer-len {params} is ignored by when --special is used
        "--special {wildcards.reference} --db {output} > {log}; " ## reference can be one of "greengenes", "silva", "rdp"
        "bracken-build -t {threads} -k {params} -l 1451"
        "-d {output} " ## median read length after filtering in 20191007_1559_MN31344_FAK76605_2bf006ff
        ">> {log}"

rule kmer_mapping_filtered:
    input:
        fastq="01_processed_data/02_trimming_filtering/{run}/{barc}/filtered.fastq", 
        target=expand("METADATA/Reference_Sequences/kraken2_{reference}", reference=config["align"]["reference"])
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
        target=expand("METADATA/Reference_Sequences/kraken2_{reference}", reference=config["align"]["reference"])
    output:
        rankS="02_analysis_results/03_kmer_mapping/{run}/{barc}/bracken/Species.bracken",
        tableS="02_analysis_results/03_kmer_mapping/{run}/{barc}/bracken/Species.kreport2",
        rankG="02_analysis_results/03_kmer_mapping/{run}/{barc}/bracken/Genus.bracken",
        tableG="02_analysis_results/03_kmer_mapping/{run}/{barc}/bracken/Genus.kreport2",
        rankF="02_analysis_results/03_kmer_mapping/{run}/{barc}/bracken/Family.bracken",
        tableF="02_analysis_results/03_kmer_mapping/{run}/{barc}/bracken/Family.kreport2"
    log:
        "02_analysis_results/03_kmer_mapping/{run}/{barc}/bracken/MeBaPiNa_ranking.log"
    benchmark:
        "02_analysis_results/03_kmer_mapping/{run}/{barc}/bracken/MeBaPiNa_ranking.benchmark.tsv"
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
        "mv 01_processed_data/03_kmer_mapping/{wildcards.run}/{wildcards.barc}/filtered_bracken.kreport2 {output.tableS}; "
        "bracken {params} "
        "-l G "
        "-d {input.target} "
        "-i {input.report} "
        "-o {output.rankG} >> {log} 2>&1; "
        "mv 01_processed_data/03_kmer_mapping/{wildcards.run}/{wildcards.barc}/filtered_bracken.kreport2 {output.tableG}; "
        "bracken {params} "
        "-l F "
        "-d {input.target} "
        "-i {input.report} "
        "-o {output.rankF} >> {log} 2>&1; "
        "mv 01_processed_data/03_kmer_mapping/{wildcards.run}/{wildcards.barc}/filtered_bracken.kreport2 {output.tableF} "

rule plot_krona_kraken2:
    input:
        "01_processed_data/03_kmer_mapping/{run}/{barc}/filtered.kraken2"
    output:
        "02_analysis_results/03_kmer_mapping/{run}/{barc}/krona/filtered.html"
    conda:
        "../envs/krona.yml"
    params:
        "-i", ## Include a wedge for queries with no hits.
        "-q 2", ## Column of input files to use as query ID. Required if magnitude files are specified. [Default: '1']
        "-t 3" ## Column of input files to use as taxonomy ID. [Default: '2']
        "-tax METADATA/Reference_Sequences/krona" ## Path to directory containing a taxonomy database to use. [Default: '/opt/miniconda/miniconda3/envs/krona/opt/krona/taxonomy']
    shell:
        "ktImportTaxonomy -i -q 2 -t 3 {input} -o {output}"

# [-s <integer>]   Column of input files to use as score. [Default: '3']
# [-m <integer>]   Column of input files to use as magnitude. If magnitude files are specified, their magnitudes will override those in this column.
