###############
## REFERENCE ##
###############

rule stat_refseq_lenstat:
    input:
        "{tmp}METADATA/Reference_Sequences/{reference}/reference.fasta.fai"
    output:
        length_stat="{tmp}03_report/Reference_Sequences/{reference}/reference_lengthdist.tsv",
        length_plot="{tmp}03_report/Reference_Sequences/{reference}/reference_lengthdist.pdf"
    conda:
        "../envs/r-diversity.yml"
    script:
        "../scripts/fai_lenstat.R"

rule stat_refseq_taxaranks:
    input:
        "{tmp}METADATA/Reference_Sequences/{reference}/krona/Species/taxlist.txt"
    output:
        "{tmp}03_report/Reference_Sequences/{reference}/reference_taxaranks.tsv"
    shell:
        "awk -F \"\\t\" '$3==\"root\"{{cnt[\"root\"]++;next}}; $3==\"domain\"{{cnt[\"domain\"]++;next}}; $3==\"phylum\"{{cnt[\"phylum\"]++;next}}; $3==\"class\"{{cnt[\"class\"]++;next}}; $3==\"order\"{{cnt[\"order\"]++;next}}; $3==\"family\"{{cnt[\"family\"]++;next}}; $3==\"genus\"{{cnt[\"genus\"]++;next}}; $3==\"species\"{{cnt[\"species\"]++;next}}; {{cnt[\"other\"]++}};"
        "END{{ printf \"%d\\troot\\n%d\\tdomain\\n%d\\tphylum\\n%d\\tclass\\n%d\\torder\\n%d\\tfamily\\n%d\\tgenus\\n%d\\tspecies\\n%d\\tother\\n\","
        "cnt[\"root\"],cnt[\"domain\"],cnt[\"phylum\"],cnt[\"class\"],cnt[\"order\"],cnt[\"family\"],cnt[\"genus\"],cnt[\"species\"],cnt[\"other\"] }}' "
        "{input} > {output}"

###############
## RAW READS ##
###############

rule stat_general_readbasecount:
    input:
        raw="{tmp}00_raw_data/{run}/report.md",
        basecall="{tmp}01_processed_data/01_basecalling/{run}/sequencing_summary.txt",
        trim="{tmp}01_processed_data/02_trimming_filtering/{run}/{barc}/trimmed.fastq.fai",
        filter="{tmp}02_analysis_results/02_trimming_filtering/{run}/nanocomp/NanoStats.txt"
    output:
        "{tmp}03_report/{timepoint}/{sample}/{run}-{barc}/read_base_counts.tsv"
    shell:
        ## get number of raw: reads
        "echo -e \"$(tail -n 4 {input.raw} | head -n 1 | cut -d\",\" -f2)\\traw_read_count_all_samples\" > {output}; "
        ## get number of basecalled: reads and bases and mean length and quality
        "awk 'BEGIN{{ all_base=0; all_qual=0 }}; "
        "NR==1{{ for(i=1; i<=NF; i++){{ " ## in first line
        "if($i == \"sequence_length_template\"){{len_col=i}} " ## get index of length column
        "else if($i == \"mean_qscore_template\"){{qual_col=i}} }}; next }}; " ## get index of quality column
        "{{ all_base = all_base + $len_col; all_qual = all_qual + $qual_col }}; " ## cumulative sum of bases and qualities
        "END{{ "
        "print (NR-1)\"\\tbasecalled_read_count_all_samples\\n\"" ## print number of reads - header
        "(all_base)\"\\tbasecalled_base_count_all_samples\\n\"" ## number of all bases
        "(all_base/(NR-1))\"\\tbasecalled_mean_length_all_samples\\n\"" ## mean bases per read
        "all_qual/(NR-1)\"\\tbasecalled_mean_quality_all_samples\"" ## mean quality
        "}}' {input.basecall} >> {output}; "
        ## get number of demultiplexed and trimmed: reads and mean read length for barcode
        "awk 'BEGIN{{ all_base = 0 }}; "
        "{{ all_base = all_base + $2 }}; " ## cumulative sum of bases
        "END{{ "
        "print (NR-1)\"\\ttrim_read_count\\n\"" ## print number of reads - header
        "(all_base)\"\\ttrim_base_count\\n\"" ## number of all bases
        "(all_base/(NR-1))\"\\ttrim_mean_length\"" ## mean bases per read
        "}}' {input.trim} >> {output}; "
        ## get number of filtered: reads and bases and mean length and quality for barcode
        "awk -v barc=\"{wildcards.barc}\" 'BEGIN{{rd_cnt=0;bp_cnt=0;rd_len=0;rd_qul=0}}; "
        "NR==1{{" ## in first row
        "for(i=1; i<=NF; i++){{" ## in all columns
        "if($i == barc){{barc_col=i}}" ## find barcode column and remember index
        "}}}}; "
        "$1==\"Number\"{{" ## in "Number of reads" row
        "gsub(\",\",\"\",$(barc_col+1)); print $(barc_col+1)\"\\tfilter_read_count\"}}; " ## remove comma and print barcode value (Note that the row name can shift the barcode index)
        "$1==\"Total\"{{" ## in "Total bases" row
        "gsub(\",\",\"\",$barc_col); print $barc_col\"\\tfilter_base_count\"}}; "
        "$1==\"Mean\"&&$3==\"length:\"{{"
        "gsub(\",\",\"\",$(barc_col+1)); print $(barc_col+1)\"\\tfilter_mean_length\"}}; "
        "$1==\"Mean\"&&$3==\"quality:\"{{"
        "gsub(\",\",\"\",$(barc_col+1)); print $(barc_col+1)\"\\tfilter_mean_quality\"}}' "
        " {input.filter} >> {output}"

#########
## OTU ##
#########

rule stat_otu_feature:
    input:
        otutable="{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/cluster_ftable.qza" ## also dummy for MeBaPiNa_q2otupick.log
    output:
        report="{tmp}03_report/{timepoint}/{sample}/{run}-{barc}/03_otu_picking-{reference}-feature_counts.tsv"
    conda:
        "../envs/qiime2.yml"
    shell:
        ## path to directory with files
        "otu_dir={input.otutable}; otu_dir=\"${{otu_dir/cluster_ftable.qza/}}\"; "
        "sleep 1; " ## make sure the log file was created
        ## get number of total: clusters and reads and mean reads per cluster (abundance)
        "if [[ ! -f \"${{otu_dir}}cluster_ftable/feature-table.biom\" ]]; then " ## if file doesnt exsist...
        "qiime tools export --input-path \"${{otu_dir}}cluster_ftable.qza\" --output-path \"${{otu_dir}}cluster_ftable\"; fi; " ## ...create it
        "awk '$2==\"observations:\"{{gsub(\",\",\"\",$3); clst_cnt=$3}}; "
        "$2==\"count:\"{{gsub(\",\",\"\",$3); read_cnt=$3}}; "
        "END {{ if(clst_cnt==0){{mean_cnt=0}}else{{mean_cnt=(read_cnt/clst_cnt)}}; "
        "printf \"%d\\ttotal_cluster_count\\n%d\\ttotal_read_count\\n%d\\ttotal_mean_abund\\n\","
        "clst_cnt,read_cnt,mean_cnt }}' "
        "<(biom summarize-table -i \"${{otu_dir}}cluster_ftable/feature-table.biom\") > {output.report}; "
        ## get number of de-novo: clusters and reads and mean reads per cluster (abundance)
        "awk 'BEGIN{{prnt_once=0}}"
        "$1==\"Clusters:\"{{clst_cnt=$2}}; "
        "$2==\"nt\"&&$3==\"in\"&&prnt_once==0{{read_cnt=$4;prnt_once++}}; "
        "END {{ if(clst_cnt==0){{mean_cnt=0}}else{{mean_cnt=(read_cnt/clst_cnt)}}; "
        "printf \"%d\\tdenovo_cluster_count\\n%d\\tdenovo_read_count\\n%d\\tdenovo_mean_abund\\n\","
        "clst_cnt,read_cnt,mean_cnt }}' "
        "\"${{otu_dir}}MeBaPiNa_q2otupick.log\" >> {output.report}; "
        ## get number of total after filtering: clusters and reads and mean reads per cluster (abundance)
        "if [[ ! -f \"${{otu_dir}}filt_ftable/feature-table.biom\" ]]; then "
        "qiime tools export --input-path \"${{otu_dir}}filt_ftable.qza\" --output-path \"${{otu_dir}}filt_ftable\"; fi; "
        "awk '$2==\"observations:\"{{gsub(\",\",\"\",$3); clst_cnt=$3}}; "
        "$2==\"count:\"{{gsub(\",\",\"\",$3); read_cnt=$3}}; "
        "END {{ if(clst_cnt==0){{mean_cnt=0}}else{{mean_cnt=(read_cnt/clst_cnt)}}; "
        "printf \"%d\\tfilter_cluster_count\\n%d\\tfilter_read_count\\n%d\\tfilter_mean_abund\\n\","
        "clst_cnt,read_cnt,mean_cnt }}' "
        "<(biom summarize-table -i \"${{otu_dir}}filt_ftable/feature-table.biom\") >> {output.report}"

rule stat_otu_taxa:
    input:
        "{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}_{reftype}/filtered.kreport2"
    output:
        report="{tmp}03_report/{timepoint}/{sample}/{run}-{barc}/03_otu_picking-{reference}_{reftype}-taxa_counts.tsv"
    shell:
        "awk 'BEGIN{{cnt_stax=0;cnt_tax=0;cnt_sfeat=0;" ## initialize variables
        "lw_rnk=substr(\"{wildcards.reftype}\",1,1)}}; " ## get first character of reference type (e.g. "S" for Species)
        "$6==\"root\"{{ cnt_feat=$2 }}; " ## get number of reads mapped to all taxa (i.e. root taxa and below)
        "$3!=0{{ cnt_tax++ }}; " ## if reads are mapped to this taxon, add one to the total number of taxa
        "$3!=0&&$4==lw_rnk{{cnt_stax++;cnt_sfeat=cnt_sfeat+$3}}; " ## ... if the current rank also is the reference type, add one here as well
        "END {{ if(cnt_tax==0){{mean_cnt=0}}else{{mean_cnt=(cnt_feat/cnt_tax)}}; " ## ensure no division by zero
        "if(cnt_stax==0){{mean_scnt=0}}else{{mean_scnt=(cnt_sfeat/cnt_stax)}}; "
        "printf \"%d\\ttotal_taxa_count\\n%d\\ttotal_read_count\\n%.2f\\ttotal_mean_abund\\n" ## print
        "%d\\t{wildcards.reftype}_taxa_count\\n%d\\t{wildcards.reftype}_read_count\\n%.2f\\t{wildcards.reftype}_mean_abund\\n\","
        "cnt_tax,cnt_feat,mean_cnt,cnt_stax,cnt_sfeat,mean_scnt }}' {input} > {output}"

rule stat_otu_taxa_diversity:
    input:
        taxlist="{tmp}METADATA/Reference_Sequences/{reference}/krona/{reftype}/taxlist.txt",
        sample_file="{tmp}02_analysis_results/03_otu_picking/{run}/{barc}/{reference}_{reftype}/kmer.counttaxlist"
    output:
        report="{tmp}03_report/{timepoint}/{sample}/{run}-{barc}/03_otu_picking-{reference}_{reftype}-taxa_diversity.tsv",
        covdist_plot="{tmp}03_report/{timepoint}/{sample}/{run}-{barc}/03_otu_picking-{reference}_{reftype}-taxa_covdist.pdf"
    conda:
        "../envs/r-diversity.yml"
    script:
        "../scripts/phyloseq_vegan.R"

###############
## ALIGNMENT ##
###############

rule stat_align_readsrates:
    input:
        "{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}/pycoqc.json"
    output:
        report="{tmp}03_report/{timepoint}/{sample}/{run}-{barc}/03_alignment-{reference}-alignment_rates.tsv"
    shell:
        "awk 'BEGIN{{prnt=0;cnt=0;strt_cnt=0}};" ## initialize variables
        "/\"alignment\"/{{prnt=1}};" ## inside the alignment section...
        "prnt==1&&/reads_number/{{gsub(\",\",\"\",$2); reads=$2; prnt=0}};" ## ... remove thousands comma and save the read number
        "/\"identity_freq_percentiles\"/{{strt_cnt=1}};" ## inside the distribution table of identity frequencies...
        "strt_cnt==1{{cnt++}};" ## ...count lines...
        "cnt==51{{gsub(\",\",\"\",$1); ident=$1}};" ## ... until median value and save it
        "/\"insertion_rate\"/{{gsub(\",\",\"\",$2); insrt=$2}};" ## save...
        "/\"deletion_rate\"/{{gsub(\",\",\"\",$2); dltn=$2}};" ## ...more...
        "/\"mismatch_rate\"/{{gsub(\",\",\"\",$2); mism=$2}};" ## ...variables
        "END{{printf \"%d\\taligned_read_count\\n%d\\tmedian_reference_identity\\n%.2f\\tdeletion_rate\\n%.2f\\tmismatch_rate\\\\n\"," ## print
        "reads,ident,insrt,dltn,mism}}' {input} > {output}"

rule stat_align_taxa:
    input:
        "{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}_{reftype}/aligned.counttaxlist"
    output:
        report="{tmp}03_report/{timepoint}/{sample}/{run}-{barc}/03_alignment-{reference}_{reftype}-taxa_counts.tsv"
    shell:
        "awk -F\"\\t\" 'BEGIN{{ cnt_stax=0;cnt_tax=0;cnt_sfeat=0;cnt_feat=0; "
        "lw_rnk=substr(\"{wildcards.reftype}\",1,1);"
        "if(lw_rnk==\"S\"){{nf_rnk=8}}else if(lw_rnk==\"G\"){{nf_rnk=7}} }}; " #!# when Species, expect 8 columns, when Genus, expect 7
        "NF==nf_rnk{{cnt_stax++;cnt_sfeat=cnt_sfeat+$1}}; " ## when number of columns matches expected number...
        "{{cnt_tax++;cnt_feat=cnt_feat+$1}}; " ## ...count to reference type
        "END {{ if(cnt_tax==0){{mean_cnt=0}}else{{mean_cnt=(cnt_feat/cnt_tax)}}; "
        "if(cnt_stax==0){{mean_scnt=0}}else{{mean_scnt=(cnt_sfeat/cnt_stax)}}; "
        "printf \"%d\\ttotal_taxa_count\\n%d\\ttotal_read_count\\n%.2f\\ttotal_mean_abund\\n"
        "%d\\t{wildcards.reftype}_taxa_count\\n%d\\t{wildcards.reftype}_read_count\\n%.2f\\t{wildcards.reftype}_mean_abund\\n\","
        "cnt_tax,cnt_feat,mean_cnt,cnt_stax,cnt_sfeat,mean_scnt }}' {input} > {output}"

rule stat_align_taxa_diversity:
    input:
        taxlist="{tmp}METADATA/Reference_Sequences/{reference}/krona/{reftype}/taxlist.txt",
        sample_file="{tmp}02_analysis_results/03_alignment/{run}/{barc}/{reference}_{reftype}/aligned.counttaxlist"
    output:
        report="{tmp}03_report/{timepoint}/{sample}/{run}-{barc}/03_alignment-{reference}_{reftype}-taxa_diversity.tsv",
        covdist_plot="{tmp}03_report/{timepoint}/{sample}/{run}-{barc}/03_alignment-{reference}_{reftype}-taxa_covdist.pdf"
    conda:
        "../envs/r-diversity.yml"
    script:
        "../scripts/phyloseq_vegan.R"

###################
## K-MER MAPPING ##
###################

rule stat_kmer_taxa:
    input:
        "{tmp}01_processed_data/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/filtered.kreport2"
    output:
        report="{tmp}03_report/{timepoint}/{sample}/{run}-{barc}/03_kmer_mapping-{reference}_{reftype}-taxa_counts.tsv"
    shell:
        "awk 'BEGIN{{cnt_stax=0;cnt_tax=0;cnt_sfeat=0;"
        "lw_rnk=substr(\"{wildcards.reftype}\",1,1)}}; "
        "$6==\"root\"{{ cnt_feat=$2 }}; "
        "$3!=0{{ cnt_tax++ }}; "
        "$3!=0&&$4==lw_rnk{{cnt_stax++;cnt_sfeat=cnt_sfeat+$3}}; "
        "END {{ if(cnt_tax==0){{mean_cnt=0}}else{{mean_cnt=(cnt_feat/cnt_tax)}}; "
        "if(cnt_stax==0){{mean_scnt=0}}else{{mean_scnt=(cnt_sfeat/cnt_stax)}}; "
        "printf \"%d\\ttotal_taxa_count\\n%d\\ttotal_read_count\\n%.2f\\ttotal_mean_abund\\n"
        "%d\\t{wildcards.reftype}_taxa_count\\n%d\\t{wildcards.reftype}_read_count\\n%.2f\\t{wildcards.reftype}_mean_abund\\n\","
        "cnt_tax,cnt_feat,mean_cnt,cnt_stax,cnt_sfeat,mean_scnt }}' {input} > {output}"

rule stat_kmer_retaxa:
    input:
        "{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/{reftype}.kreport2"
    output:
        report="{tmp}03_report/{timepoint}/{sample}/{run}-{barc}/03_kmer_mapping-{reference}_{reftype}-retaxa_counts.tsv"
    shell:
        "awk 'BEGIN{{cnt_stax=0;cnt_tax=0;cnt_sfeat=0;"
        "lw_rnk=substr(\"{wildcards.reftype}\",1,1)}}; "
        "$6==\"root\"{{ cnt_feat=$2 }}; "
        "$3!=0{{ cnt_tax++ }}; "
        "$3!=0&&$4==lw_rnk{{cnt_stax++;cnt_sfeat=cnt_sfeat+$3}}; "
        "END {{ if(cnt_tax==0){{mean_cnt=0}}else{{mean_cnt=(cnt_feat/cnt_tax)}}; "
        "if(cnt_stax==0){{mean_scnt=0}}else{{mean_scnt=(cnt_sfeat/cnt_stax)}}; "
        "printf \"%d\\ttotal_taxa_count\\n%d\\ttotal_read_count\\n%.2f\\ttotal_mean_abund\\n"
        "%d\\t{wildcards.reftype}_taxa_count\\n%d\\t{wildcards.reftype}_read_count\\n%.2f\\t{wildcards.reftype}_mean_abund\\n\","
        "cnt_tax,cnt_feat,mean_cnt,cnt_stax,cnt_sfeat,mean_scnt }}' {input} > {output}"

rule stat_kmer_taxa_diversity:
    input:
        taxlist="{tmp}METADATA/Reference_Sequences/{reference}/krona/{reftype}/taxlist.txt",
        sample_file="{tmp}02_analysis_results/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/kmer.counttaxlist"
    output:
        report="{tmp}03_report/{timepoint}/{sample}/{run}-{barc}/03_kmer_mapping-{reference}_{reftype}-taxa_diversity.tsv",
        covdist_plot="{tmp}03_report/{timepoint}/{sample}/{run}-{barc}/03_kmer_mapping-{reference}_{reftype}-taxa_covdist.pdf"
    conda:
        "../envs/r-diversity.yml"
    script:
        "../scripts/phyloseq_vegan.R"
