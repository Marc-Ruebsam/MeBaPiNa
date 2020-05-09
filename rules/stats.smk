###############
## REFERENCE ##
###############

rule stat_refseq_lendist:
    input:
        "{tmp}METADATA/Reference_Sequences/{reference}/reference.fasta"
    output:
        "{tmp}.tmp/stats/{reference}/refseq_lendist.tsv"
    conda:
        "../envs/samtools.yml"
    shell:
        "samtools faidx {input}; "
        "cut -f 2 {input}.fai | sort -n | uniq -c | awk '{{print $1\"\t\"$2}}' > {output}"

sort -nk2 reference.fasta.fai | awk 'BEGIN{c=0};
    {cnt_sort[c++]=$2};
    END{
        sm = 0;
        for(key in cnt_sort){ sm = sm + cnt_sort[key] };
        men = sm / c;
        if((c % 2) == 1){ medn = cnt_sort[ int(c/2) ];
        }else{ medn = ( cnt_sort[c/2] + cnt_sort[c/2-1] ) / 2 };
        print medn"\t"men }'

## taxonomic ranks and frequency
# cut -f 3 taxlist.txt | sort | uniq -c | sort -nr | egrep " root| domain| phylum| class| order| family| genus| species";
# cut -f 3 taxlist.txt | sort | uniq -c | sort -nr | egrep -v " root| domain| phylum| class| order| family| genus| species" | awk 'BEGIN{cnt=0};{cnt=cnt+$1};END{print cnt" other"}'

###############
## RAW READS ##
###############

rule stat_raw_read_count:
    input:
        "{tmp}00_raw_data/{run}/report.md"
    output:
        "{tmp}.tmp/stats/{run}/raw_read_count.tsv"
    shell:
        "echo -e \"$(tail -n 4 {input} | head -n 1 | cut -d\",\" -f2)\\traw_read_count\" >> {output}"

###################
## DEMULTIPLEXED ##
###################

rule stat_demux_assigned_counts:
    input:
        "{tmp}02_analysis_results/01_basecalling/{run}/nanocomp/NanoStats.txt"
    output:
        "{tmp}.tmp/stats/{run}/demux_assigned_counts.tsv"
    shell:
        "awk 'BEGIN{{rd_cnt=0;bp_cnt=0;rd_len=0;rd_qul=0}}; "
        "NR==1{{"
        "for(i=1; i<=NF; i++){{"
        "if($i == \"unclassified\"){{uncl_col=i}}"
        "}}}}; "
        "$1==\"Number\"{{"
        "for(i=4; i<NF; i++){{"
        "gsub(\",\",\"\",$i); rd_cnt = rd_cnt + $i"
        "}}}}; "
        "$1==\"Total\"{{"
        "for(i=3; i<NF; i++){{"
        "gsub(\",\",\"\",$i); bp_cnt = bp_cnt + $i"
        "}}}}; "
        "$1==\"Median\"&&$3==\"length:\"{{"
        "for(i=4; i<NF; i++){{"
        "gsub(\",\",\"\",$i); rd_len = rd_len + $i"
        "}}; rd_len = rd_len / (NF-3)}}; "
        "$1==\"Median\"&&$3==\"quality:\"{{"
        "for(i=4; i<NF; i++){{"
        "gsub(\",\",\"\",$i); rd_qul = rd_qul + $i"
        "}}; rd_qul = rd_qul / (NF-3)}}; "
        "END{{print "
        "rd_cnt\"\\tassigned_read_count\\n\""
        "bp_cnt\"\\tassigned_base_count\\n\""
        "rd_len\"\\tassigned_mean_length\\n\""
        "rd_qul\"\\tassigned_mean_quality\""
        "}}' {input} >> {output}"

# filtered_reads <- c(
#   44,
#   30,
#   673100,
#   1487303,
#   718362,
#   857069,
#   678881,
#   530264,
#   1341020,
#   988081,
#   724721,
#   625320,
#   66260,
#   283913,
#   723653,
#   191520,
#   192578,
#   157178,
#   138630,
#   61976,
#   78003,
#   158098,
#   230634,
#   932741
# )

###########
## K-MER ##
###########

## run time
# find -name "MeBaPiNa_kmermap_filtered.benchmark.tsv" -exec tail -n 1 {} \; | sort -nr | head -n 1
## number of total and species taxa, total ans species reads and median reads per Species
for fl in $(find 16S_Metabarcoding/01_processed_data/03_kmer_mapping/ -name "filtered.kreport2" | sort)
    do
    awk '
    BEGIN{
        c=0;
        cnt_stax=0;
        cnt_htax=0 };
    $6=="root"{ cnt_hfeat=$2 };
    $3!=0{ cnt_htax++ };
    $3!=0&&$4=="S"{
        cnt_stax++;
        cnt_sfeat[c++]=$3 };
    END{
        sm = 0;
        for(key in cnt_sfeat){ sm = sm + cnt_sfeat[key] };
        if((c % 2) == 1){ medn = cnt_sfeat[ int(c/2) ];
        }else{ medn = ( cnt_sfeat[c/2] + cnt_sfeat[c/2-1] ) / 2 };
        print cnt_htax"\t"cnt_stax"\t"cnt_hfeat"\t"sm"\t"medn"\t" }
    ' <(sort -nk3 $fl)
done
## number of species taxa, species reads and median reads per Species
for fl in $(find 16S_Metabarcoding/02_analysis_results/03_kmer_mapping/ -name "Species.kreport2" | sort)
    do
    awk '
    BEGIN{
        c=0;
        cnt_stax=0;
        cnt_htax=0 };
    $6=="root"{ cnt_hfeat=$2 };
    $3!=0{ cnt_htax++ };
    $3!=0&&$4=="S"{
        cnt_stax++;
        cnt_sfeat[c++]=$3 };
    END{
        sm = 0;
        for(key in cnt_sfeat){ sm = sm + cnt_sfeat[key] };
        if((c % 2) == 1){ medn = cnt_sfeat[ int(c/2) ];
        }else{ medn = ( cnt_sfeat[c/2] + cnt_sfeat[c/2-1] ) / 2 };
        print cnt_stax"\t"sm"\t"medn"\t" }
        # print cnt_htax"\t"cnt_stax"\t"cnt_hfeat"\t"sm"\t"medn"\t" } #!!!# taxa are the same, reads not
    ' <(sort -nk3 $fl)
done

###############
## ALIGNMENT ##
###############

rule stat_ref_coverage:
    input:
        bam=expand("{tmp}01_processed_data/03_alignment/{run}/{barc}/{reference}/{altype}_filteredsorted.bam",
        tmp = config["experiments"]["tmp"], run = RUNS, barc = SAMPLES.keys(), reference = config['reference']['source'], reftype = config['reference']['rank'], altype = "aligned"),
        reference="{tmp}METADATA/Reference_Sequences/{reference}/reference.fasta"
    output:
        "{tmp}.tmp/stats/{run}/{reference}/refseq_coverage.tsv"
    conda:
        "../envs/samtools.yml"
    params:
        "-a",
        "-l " + config['filtering']['len_min'] ## should be redundant
    shell:
        "samtools depth {params} "
        "--reference {input.referemce} {input.bam} "
        "> {output}"

## time
# find -name "MeBaPiNa_aligning_filtered.benchmark.tsv" -exec tail -n 1 {} \; | sort -nr | head -n 1
## aligned reads, insertion, deletion and mismatch rate and read identity frequency
for fl in $(find 16S_Metabarcoding/02_analysis_results/03_alignment/ -name "pycoqc.json" | sort);
    do
    awk '
    BEGIN{prnt=0;cnt=0;strt_cnt=0};
    /"alignment"/{prnt=1};
    /identity_freq_percentiles/{strt_cnt=1};
    strt_cnt==1{cnt++};
    cnt==51{gsub(",","",$1); ident=$1};
    prnt==1&&/reads_number/{gsub(",","",$2); reads=$2; prnt=0};
    /"insertion_rate"/{gsub(",","",$2); insrt=$2};
    /"deletion_rate"/{gsub(",","",$2); dltn=$2};
    /"mismatch_rate"/{gsub(",","",$2); mism=$2};
    END{print reads"\t"insrt"\t"dltn"\t"mism"\t"ident}
    ' $fl;
done
## faster then samtools view -c
# samtools flagstat
# samtools stats
## total taxa (all species), assiged reads, median reads per taxa
for fl in $(find 16S_Metabarcoding/02_analysis_results/03_alignment/ -name "aligned.counttaxlist" | sort)
    do
    awk '
    BEGIN{c=0};
    {cnt_sfeat[c++]=$1};
    END{
        sm = 0;
        for(key in cnt_sfeat){ sm = sm + cnt_sfeat[key] };
        if((c % 2) == 1){ medn = cnt_sfeat[ int(c/2) ];
        }else{ medn = ( cnt_sfeat[c/2] + cnt_sfeat[c/2-1] ) / 2 };
        print c"\t"sm"\t"medn"\t" }
    ' <(sort -nk1 $fl)
done

#########
## OTU ##
#########

rule stat_otu:
    input:
        otutable="{tmp}01_processed_data/03_otu_picking/{run}/{barc}/{reference}/cluster_ftable.qza" ## also dummy for MeBaPiNa.log
    output:
        "{tmp}.tmp/stats/{run}/{barc}/{reference}/otu_counts_counts_{barc}.tsv"
    shell:
        "otu­_dir={input.otutable}; otu_dir=\"${{otu_dir/cluster_ftable.qza/}}\" > {log} 2>&1; "
        "sleep 1;" ## make sure the log file was created
        "awk '$1==\"Clusters:\"{{print \"denovo_cluster\\t\"$2}}' \"${{otu_dir}}MeBaPiNa_q2otupick.log\" >> {output}; "
        "qiime tools export -i \"${{otu_dir}}cluster_ftable.qza\" -o \"${{otu_dir}}cluster_ftable\" > {log} 2>&1; "
        "awk '$2==\"observations:\"{{print \"total otu clusters\\t\"$3}}; $2==\"count:\"{{print \"total otu reads\\t\"$3}}'
        "<(biom summarize-table -i \"${{otu_dir}}cluster_ftable\feature-table.biom\"); done"

## time
find -name "MeBaPiNa_kmermap_q2rereplicate.benchmark.tsv" -exec tail -n 1 {} \; | sort -nr | head -n 1
## de novo clusters
for fl in $(find 16S_Metabarcoding/01_processed_data/03_otu_picking/ -name "MeBaPiNa_q2otupick.log" | sort)
    do
    awk '$1=="Clusters:"{print $2}' $fl
done

## number of features, numbr of reads and median reads per feature BEFORE filterng
for fl in $(find 16S_Metabarcoding/01_processed_data/03_otu_picking/ -name "cluster_ftable.qza" | sort)
    do
    qiime tools export --input-path ${fl} --output-path "${fl/.qza/}" > /dev/null
    biom convert -i "${fl/.qza/}/feature-table.biom" -o "${fl/.qza/}/feature-table.tsv" --to-tsv > /dev/null
    awk '
    BEGIN{c=0};
    !/^\#/{cnt_sfeat[c++]=$2};
    END{
        sm = 0;
        for(key in cnt_sfeat){ sm = sm + cnt_sfeat[key] };
        if((c % 2) == 1){ medn = cnt_sfeat[ int(c/2) ];
        }else{ medn = ( cnt_sfeat[c/2] + cnt_sfeat[c/2-1] ) / 2 };
        print c"\t"sm"\t"medn"\t" }
    ' <(sort -nk2 "${fl/.qza/}/feature-table.tsv")
done

## number of features, numbr of reads and median reads per feature AFTER filterng
for fl in $(find 16S_Metabarcoding/01_processed_data/03_otu_picking/ -name "filt_ftable.qza" | sort)
    do
    qiime tools export --input-path ${fl} --output-path "${fl/.qza/}" > /dev/null
    biom convert -i "${fl/.qza/}/feature-table.biom" -o "${fl/.qza/}/feature-table.tsv" --to-tsv > /dev/null
    awk '
    BEGIN{c=0};
    !/^\#/{cnt_sfeat[c++]=$2};
    END{
        sm = 0;
        for(key in cnt_sfeat){ sm = sm + cnt_sfeat[key] };
        if((c % 2) == 1){ medn = cnt_sfeat[ int(c/2) ];
        }else{ medn = ( cnt_sfeat[c/2] + cnt_sfeat[c/2-1] ) / 2 };
        print c"\t"sm"\t"medn"\t" }
    ' <(sort -nk2 "${fl/.qza/}/feature-table.tsv")
done

## time kmer
# find -name "MeBaPiNa_kmermap_q2rereplicate.benchmark.tsv" -exec tail -n 1 {} \; | sort -nr | head -n 1
for fl in $(find 16S_Metabarcoding/01_processed_data/03_otu_picking/ -name "filtered.kreport2" | sort)
    do
    awk '
    BEGIN{
        c=0;
        cnt_stax=0;
        cnt_htax=0 };
    $6=="root"{ cnt_hfeat=$2 };
    $3!=0{ cnt_htax++ };
    $3!=0&&$4=="S"{
        cnt_stax++;
        cnt_sfeat[c++]=$3 };
    END{
        sm = 0;
        for(key in cnt_sfeat){ sm = sm + cnt_sfeat[key] };
        if((c % 2) == 1){ medn = cnt_sfeat[ int(c/2) ];
        }else{ medn = ( cnt_sfeat[c/2] + cnt_sfeat[c/2-1] ) / 2 };
        print cnt_htax"\t"cnt_stax"\t"cnt_hfeat"\t"sm"\t"medn"\t" }
    ' <(sort -nk3 $fl)
done



# for fl in $(find -name "*.qza")
#   do
#   echo ${fl}
#   qiime tools export --input-path ${fl} --output-path "${fl/.qza/}"
# done
# for fl in $(find -name "*.biom")
#   do
#   echo ${fl}
#   biom summarize-table -i ${fl}
#   biom convert -i ${fl} -o ${fl/.biom/.tsv} --to-tsv
# done

# for fl in $(find -name "*.biom")
#   do
#   echo ${fl}
#   biom summarize-table -i ${fl}
#   biom convert -i ${fl} -o ${fl/.biom/.tsv} --to-tsv
# done

###################################################################################################################

####################
## mock community ##
####################

## species -> genus -> family
taxids_S=(80813 102212 118750 98728 82239 94774 89641 78858)
for id in "${taxids_S[@]}"; do awk -v taxID="$id" '$5==taxID{print $2}' Species.kreport2; done
taxids_G=(3723 1831 46463 46474 45258 45295 1688 1837)
for id in "${taxids_G[@]}"; do awk -v taxID="$id" '$5==taxID{print $2}' Species.kreport2; done
taxids_F=(3718 1828 46454 45256 45287 1680 1836)
for id in "${taxids_F[@]}"; do awk -v taxID="$id" '$5==taxID{print $2}' Species.kreport2; done


names_S=("Bacteria\tProteobacteria\tGammaproteobacteria\tPseudomonadales\tPseudomonadaceae\tPseudomonas\tPseudomonas aeruginosa" "Bacteria\tFirmicutes\tBacilli\tLactobacillales\tEnterococcaceae\tEnterococcus\tEnterococcus faecalis" "Bacteria\tProteobacteria\tGammaproteobacteria\tEnterobacterales\tEnterobacteriaceae\tEscherichia-Shigella\tEscherichia coli" "Bacteria\tProteobacteria\tGammaproteobacteria\tEnterobacterales\tEnterobacteriaceae\tSalmonella\tSalmonella enterica" "Bacteria\tFirmicutes\tBacilli\tLactobacillales\tListeriaceae\tListeria\tListeria monocytogenes" "Bacteria\tFirmicutes\tBacilli\tStaphylococcales\tStaphylococcaceae\tStaphylococcus\tStaphylococcus aureus" "Bacteria\tFirmicutes\tBacilli\tBacillales\tBacillaceae\tBacillus\tBacillus subtilis" "Bacteria\tFirmicutes\tBacilli\tLactobacillales\tLactobacillaceae\tLactobacillus\tLactobacillus fermentum")
for nm in "${names_S[@]}"; do awk -F"\t" 'BEGIN{sm=0}; {sm=sm+$1}; END{print sm}' <(grep -P "$nm" aligned.counttaxlist); done
names_G=("Bacteria\tProteobacteria\tGammaproteobacteria\tPseudomonadales\tPseudomonadaceae\tPseudomonas" "Bacteria\tFirmicutes\tBacilli\tLactobacillales\tEnterococcaceae\tEnterococcus" "Bacteria\tProteobacteria\tGammaproteobacteria\tEnterobacterales\tEnterobacteriaceae\tEscherichia-Shigella" "Bacteria\tProteobacteria\tGammaproteobacteria\tEnterobacterales\tEnterobacteriaceae\tSalmonella" "Bacteria\tFirmicutes\tBacilli\tLactobacillales\tListeriaceae\tListeria" "Bacteria\tFirmicutes\tBacilli\tStaphylococcales\tStaphylococcaceae\tStaphylococcus" "Bacteria\tFirmicutes\tBacilli\tBacillales\tBacillaceae\tBacillus" "Bacteria\tFirmicutes\tBacilli\tLactobacillales\tLactobacillaceae\tLactobacillus")
for nm in "${names_G[@]}"; do awk -F"\t" 'BEGIN{sm=0}; {sm=sm+$1}; END{print sm}' <(grep -P "$nm" aligned.counttaxlist); done
names_F=("Bacteria\tProteobacteria\tGammaproteobacteria\tPseudomonadales\tPseudomonadaceae" "Bacteria\tFirmicutes\tBacilli\tLactobacillales\tEnterococcaceae" "Bacteria\tProteobacteria\tGammaproteobacteria\tEnterobacterales" "Bacteria\tFirmicutes\tBacilli\tLactobacillales\tListeriaceae" "Bacteria\tFirmicutes\tBacilli\tStaphylococcales\tStaphylococcaceae" "Bacteria\tFirmicutes\tBacilli\tBacillales\tBacillaceae" "Bacteria\tFirmicutes\tBacilli\tLactobacillales\tLactobacillaceae")
for nm in "${names_F[@]}"; do awk -F"\t" 'BEGIN{sm=0}; {sm=sm+$1}; END{print sm}' <(grep -P "$nm" aligned.counttaxlist); done

## species -> genus -> family
taxids_S=(80813 102212 118750 98728 82239 94774 89641 78858)
for id in "${taxids_S[@]}"; do awk -v taxID="$id" '$5==taxID{print $2}' filtered.kreport2; done
taxids_G=(3723 1831 46463 46474 45258 45295 1688 1837)
for id in "${taxids_G[@]}"; do awk -v taxID="$id" '$5==taxID{print $2}' filtered.kreport2; done
taxids_F=(3718 1828 46454 45256 45287 1680 1836)
for id in "${taxids_F[@]}"; do awk -v taxID="$id" '$5==taxID{print $2}' filtered.kreport2; done









#
# rule deseq2:
#     input:
#         expand("{tmp}01_processed_data/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/Species.counttaxlist",
#         tmp = config["experiments"]["tmp"], run = RUNS, barc = SAMPLES.keys(), reference = config["reference"]["source"], reftype = config["reference"]["type"])
#     output:
#         "doesntexist"
#     shell:
#         "echo {input}"
#
#
#













# > ls()
# [1] "count_tab"       "sample_info_tab" "tax_tab"
#
# > str(count_tab)
# 'data.frame':   2498 obs. of  16 variables:
#  $ BW1  : int  0 0 0 0 0 0 0 0 0 0 ...
#  $ BW2  : int  0 0 1 0 0 0 0 1 0 0 ...
#  $ R10  : int  0 459 44 392 0 76 160 11 11 5 ...
#  $ R11BF: int  1762 1 0 0 0 4 0 3 3 0 ...
#  $ R11  : int  4 130 26 236 0 82 48 0 8 4 ...
#  $ R12  : int  4 78 148 6 0 74 117 77 76 71 ...
#  $ R1A  : int  0 43 204 2 0 64 71 115 183 70 ...
#  $ R1B  : int  0 38 153 0 0 74 81 80 258 102 ...
#  $ R2   : int  3 90 154 38 0 115 125 104 72 118 ...
#  $ R3   : int  0 148 68 64 0 128 90 94 41 116 ...
#  $ R4   : int  5 159 113 62 0 155 166 213 111 284 ...
#  $ R5   : int  7 71 140 56 0 108 78 107 150 154 ...
#  $ R6   : int  4 8 193 0 0 122 29 77 129 87 ...
#  $ R7   : int  2 0 5 0 0 81 0 0 49 48 ...
#  $ R8   : int  2 178 39 260 0 38 123 98 0 0 ...
#  $ R9   : int  0 159 19 160 1199 56 29 134 11 5 ...
#  > head(count_tab)
#        BW1 BW2 R10 R11BF R11 R12 R1A R1B  R2  R3  R4  R5  R6 R7  R8   R9
#  ASV_1   0   0   0  1762   4   4   0   0   3   0   5   7   4  2   2    0
#  ASV_2   0   0 459     1 130  78  43  38  90 148 159  71   8  0 178  159
#  ASV_3   0   1  44     0  26 148 204 153 154  68 113 140 193  5  39   19
#  ASV_4   0   0 392     0 236   6   2   0  38  64  62  56   0  0 260  160
#  ASV_5   0   0   0     0   0   0   0   0   0   0   0   0   0  0   0 1199
#  ASV_6   0   0  76     4  82  74  64  74 115 128 155 108 122 81  38   56
#
# > str(tax_tab)
#  chr [1:2498, 1:6] "Bacteria" "Bacteria" "Bacteria" "Bacteria" "Bacteria" ...
#  - attr(*, "dimnames")=List of 2
#   ..$ : chr [1:2498] "ASV_1" "ASV_2" "ASV_3" "ASV_4" ...
#   ..$ : chr [1:6] "Kingdom" "Phylum" "Class" "Order" ...
#  > head(tax_tab)
#        Kingdom    Phylum           Class                 Order
#  ASV_1 "Bacteria" "Proteobacteria" "Alphaproteobacteria" "Kordiimonadales"
#  ASV_2 "Bacteria" "Nitrospirae"    "Nitrospira"          "Nitrospirales"
#  ASV_3 "Bacteria" "Proteobacteria" "Alphaproteobacteria" "Rhodovibrionales"
#  ASV_4 "Bacteria" "Proteobacteria" "Gammaproteobacteria" "pItb-vmat-80"
#  ASV_5 "Bacteria" "Proteobacteria" "Gammaproteobacteria" "UBA10353_marine_group"
#  ASV_6 "Bacteria" "Proteobacteria" "Gammaproteobacteria" "Steroidobacterales"
#        Family           Genus
#  ASV_1 NA               NA
#  ASV_2 "Nitrospiraceae" "Nitrospira"
#  ASV_3 "Kiloniellaceae" NA
#  ASV_4 NA               NA
#  ASV_5 NA               NA
#  ASV_6 "Woeseiaceae"    "Woeseia"
#
# > str(sample_info_tab)
# 'data.frame':   16 obs. of  4 variables:
#  $ temp : num  2 2 13.7 7.3 7.3 NA 8.6 8.6 8.6 12.7 ...
#  $ type : Factor w/ 3 levels "biofilm","rock",..: 3 3 2 1 2 2 2 2 2 2 ...
#  $ char : Factor w/ 5 levels "altered","biofilm",..: 5 5 4 2 4 1 1 1 1 1 ...
#  $ color: chr  "blue" "blue" "black" "darkgreen" ...
