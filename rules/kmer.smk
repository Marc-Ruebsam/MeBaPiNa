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
    shell:
        "out1={output[0]}; "
        "cp ${{out1/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} ${{out1}}; "
        "out2={output[1]}; "
        "cp ${{out2/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} ${{out2}}; "
        "lg={log}; "
        "cat ${{lg/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} > {log}"

rule retax_kmermap:
    input:
        report="{tmp}01_processed_data/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/filtered.kreport2",
        output="{tmp}01_processed_data/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/filtered.kraken2",
        krakdb="{tmp}METADATA/Reference_Sequences/{reference}/kraken2/{reftype}/database.kraken"
    output:
        rank="{tmp}01_processed_data/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/{reftype}.bracken",
        table="{tmp}01_processed_data/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/{reftype}.kreport2"
    log:
        "{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_retax_kmermap.log"
    benchmark:
        "{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_retax_kmermap.benchmark.tsv"
    shell:
        "out1={output[0]}; "
        "cp ${{out1/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} ${{out1}}; "
        "out2={output[1]}; "
        "cp ${{out2/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} ${{out2}}; "
        "lg={log}; "
        "cat ${{lg/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} > {log}"

rule counttax_kmermap:
    input:
        # kreport="{tmp}01_processed_data/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/filtered.kreport2", ## kraken2
        kreport="{tmp}01_processed_data/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/{reftype}.kreport2", ## bracken
        kronataxlist="{tmp}METADATA/Reference_Sequences/{reference}/krona/{reftype}/taxlist.txt"
    output:
        counttaxlist="{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/kmer.counttaxlist"
    log:
        "{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_counttax_kmermap.log"
    benchmark:
        "{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/MeBaPiNa_counttax_kmermap.benchmark.tsv"
    shell:
        "out1={output[0]}; "
        "cp ${{out1/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} ${{out1}}; "
        "lg={log}; "
        "cat ${{lg/16S_Metabarcoding/\"16S_Metabarcoding/_Temp\"}} > {log}"
