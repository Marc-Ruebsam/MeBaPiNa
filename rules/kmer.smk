## KRAKEN2 ##

rule kmermap_filtered:
    input:
        fastq="{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/filtered.fastq",
        krakdb="{tmp}METADATA/Reference_Sequences/{reference}/kraken2/{reftype}/database.kraken"
    output:
        report="{tmp}01_processed_data/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/filtered.kreport2",
        output="{tmp}01_processed_data/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/filtered.kraken2"
    log:
        "{tmp}01_processed_data/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_kmermap_filtered.log"
    benchmark:
        "{tmp}01_processed_data/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_kmermap_filtered.benchmark.tsv"
    conda:
        "../envs/kraken2.yml"
    threads:
        8
    params:
        "--confidence " + config["filtering"]["min_confidence"] ## how many of the k-mers have to map to a reference to be assigned (higher taxonomies accumulate the counts of lower ones)
    shell:
        "target={input.krakdb}; target=\"${{target/database.kraken/}}\" > {log} 2>&1; "
        "kraken2 --threads {threads} {params} "
        "--db ${{target}} "
        "--output {output.output} " ## information per sequence
        "--report {output.report} " ## information per taxon
        "{input.fastq} >> {log} 2>&1"

rule retax_kmermap:
    input:
        report="{tmp}01_processed_data/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/filtered.kreport2",
        output="{tmp}01_processed_data/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/filtered.kraken2",
        krakdb="{tmp}METADATA/Reference_Sequences/{reference}/kraken2/{reftype}/database.kraken"
    output:
        # rankF="{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/Family.bracken",
        # tableF="{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/Family.kreport2",
        rankG="{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/Genus.bracken",
        tableG="{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/Genus.kreport2",
        rankS="{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/Species.bracken",
        tableS="{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/Species.kreport2"
    log:
        "{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_retax_kmermap.log"
    benchmark:
        "{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_retax_kmermap.benchmark.tsv"
    conda:
        "../envs/kraken2.yml"
    params:
        "-t " + config['filtering']['min_featurereads'], ## [Default = 10]:: specifies the minimum number of reads required for a classification at the specified rank. Any classifications with less than the specified threshold will not receive additional reads from higher taxonomy levels when distributing reads for abundance estimation.
        "-r 1451" ## median read length after filtering in 20191007_1559_MN31344_FAK76605_2bf006ff
    shell:
        "target={input.krakdb}; target=\"${{target/database.kraken/}}\" > {log} 2>&1; "
        "in_dir={input.report}; in_dir=\"${{in_dir/filtered.kreport2/}}\" > {log} 2>&1; "
        # "bracken {params} "
        # "-l F " ## [Default = 'S', Options = 'D','P','C','O','F','G','S']:: specifies the taxonomic rank to analyze. Each classification at this specified rank will receive an estimated number of reads belonging to that rank after abundance estimation.
        # "-d ${{target}} "
        # "-i {input.report} "
        # "-o {output.rankF} >> {log} 2>&1; "
        # "mv ${{in_dir}}filtered_bracken.kreport2 {output.tableF} >> {log} 2>&1; "
        "bracken {params} "
        "-l G "
        "-d ${{target}} "
        "-i {input.report} "
        "-o {output.rankG} >> {log} 2>&1; "
        "mv ${{in_dir}}filtered_bracken.kreport2 {output.tableG} >> {log} 2>&1; "
        "bracken {params} "
        "-l S "
        "-d ${{target}} "
        "-i {input.report} "
        "-o {output.rankS} >> {log} 2>&1; "
        "mv ${{in_dir}}/filtered_bracken.kreport2 {output.tableS} >> {log} 2>&1"

rule counttax_kmermap:
    input:
        # kreport="{tmp}01_processed_data/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/filtered.kreport2", ## kraken2
        kreport="{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/{reftype}.kreport2", ## bracken
        kronataxlist="{tmp}METADATA/Reference_Sequences/{reference}/krona/{reftype}/taxlist.txt"
    output:
        counttaxlist="{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/kmer.counttaxlist"
    log:
        "{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_counttax_kmermap.log"
    benchmark:
        "{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_counttax_kmermap.benchmark.tsv"
    conda:
        "../envs/python.yml"
    script:
        "../scripts/convert_kreport.py"
