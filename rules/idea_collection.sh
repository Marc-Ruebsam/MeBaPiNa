# rule porechop:
#     input:
#         "01_processeddata/{run}/filter/{barc}.fastq"
#     output:
#         "01_processeddata/{run}/trim/{barc}.fastq"
#     log:
#         "01_processeddata/{run}/trim/{barc}_MeBaPiNa_porechop.log"
#     benchmark:
#         "01_processeddata/{run}/trim/{barc}_MeBaPiNa_porechop.benchmark.tsv"
#     conda:
#         "../envs/porechop.yml"
#     threads:
#         config["machine"]["cpu"]
#     params:
#         "--format fastq",
#         "--check_reads 1000", ## This many reads will be aligned to all possible adapters to determine which adapter sets are present (default: 10000)
#         "--barcode_threshold 75.0", ## A read must have at least this percent identity to a barcode to be binned (default: 75.0)
#         "--barcode_diff 5.0", ## If the difference between a read's best barcode identity and its second-best barcode identity is less than this value, it will not be put in a barcode bin (to exclude cases which are too close to call) (default: 5.0)
#         # "--require_two_barcodes", ## Reads will only be put in barcode bins if they have a strong match for the barcode on both their start and end (default: a read can be binned with a match at its start or end)
#         "--adapter_threshold 90.0", ## An adapter set has to have at least this percent identity to be labelled as present and trimmed off (0 to 100) (default: 90.0)
#         "--end_threshold 75.0", ## Adapters at the ends of reads must have at least this percent identity to be removed (0 to 100) (default: 75.0)
#         "--extra_end_trim 2", ## This many additional bases will be removed next to adapters found at the ends of reads (default: 2)
#         "--verbosity 1" ##  Level of progress information: 0 = none, 1 = some, 2 = lots, 3 = full - output will go to stdout if reads are saved to a file and stderr if reads are printed to stdout (default: 1)
#     shell:
#         "mkdir -p 01_processeddata/{wildcards.run}/trim/{wildcards.barc}; "
#         "porechop --threads {threads} {params} --input {input} --barcode_dir 01_processeddata/{wildcards.run}/trim/{wildcards.barc} > {log} 2>&1; "
#         "barcode_string={wildcards.barc}; "
#         "gzip --stdout 01_processeddata/{wildcards.run}/trim/{wildcards.barc}/BC${{barcode_string: -2}}.fastq > 01_processeddata/{wildcards.run}/trim/{wildcards.barc}.fastq.gz"

# samtools:
# plot-bamstats -s /ag-halama/Microbiome/16S_Metabarcoding/00_rawdata/reference_sequences/silva_arb/SILVA_132_SSURef_Nr99_tax_silva.fasta > /ag-halama/Microbiome/16S_Metabarcoding/00_rawdata/reference_sequences/silva_arb/SILVA_132_SSURef_Nr99_tax_silva.fasta.gc
# plot-bamstats -r /ag-halama/Microbiome/16S_Metabarcoding/00_rawdata/reference_sequences/silva_arb/SILVA_132_SSURef_Nr99_tax_silva.fasta.gc -p /ag-halama/Microbiome/16S_Metabarcoding/02_analysis/dummy_pool/align/samtools/ <(samtools stats unclassified_alignment.sam)
# plot-bamstats -s /ag-halama/Microbiome/16S_Metabarcoding/00_rawdata/reference_sequences/lambda/lambda_full_NC_001416.1_13-AUG-2018.fasta > /ag-halama/Microbiome/16S_Metabarcoding/00_rawdata/reference_sequences/lambda/lambda_full_NC_001416.1_13-AUG-2018.fasta.gc
# plot-bamstats -r /ag-halama/Microbiome/16S_Metabarcoding/00_rawdata/reference_sequences/lambda/lambda_full_NC_001416.1_13-AUG-2018.fasta.gc -p /ag-halama/Microbiome/16S_Metabarcoding/02_analysis/20191001_HD-T980_RapidHMW/Lamda/Lambda/20191001_1411_MN31344_FAK77797_2a466f85/align/samtools/ <(samtools stats unclassified_alignment.sam)
# samtools depth -a -l 1000 --reference /ag-halama/Microbiome/16S_Metabarcoding/00_rawdata/reference_sequences/silva_arb/SILVA_132_SSURef_Nr99_tax_silva.fasta  barcode04_alignment_sorted.bam barcode05_alignment_sorted.bam barcode06_alignment_sorted.bam > /ag-halama/Microbiome/16S_Metabarcoding/02_analysis/16S_nanopore_primer/Testpool/20191007_1559_MN31344_FAK76605_2bf006ff/align/samtools_depth.txt
# samtools flagstat

# qualimap:
# qualimap bamqc -bam barcode10_alignment_sorted.bam -hm 3 -nt 2 --java-mem-size=2G --output-genome-coverage testpicture -outdir /ag-halama/Microbiome/16S_Metabarcoding/02_analysis/16S_nanopore_primer/Testpool/20191007_1559_MN31344_FAK76605_2bf006ff/align/qualimap/
# qualimap rnaseq
# qualimap multi-bamqc
# awk 'BEGIN{prnt_flag=1}; /fastq/{next}; /Chromosome stats/{prnt_flag=0}; /summary section/{prnt_flag=1}; prnt_flag==1' qualimapReport.html > test.html

# wub:
# "bam_accuracy.py -g bam_accuracy.tsv -l bam_accuracy_reads.tsv -r bam_accuracy.pdf -e ../../../01_processeddata/dummy_pool/align_calibration_strands/lambda_alignment_sorted.bam"
# "bam_alignment_length.py -t bam_alignment_length.tsv ../../../01_processeddata/dummy_pool/align_calibration_strands/lambda_alignment_sorted.bam"
# "bam_alignment_qc.py -f ../../../00_rawdata/reference_sequences/lambda/lambda_3.6kb.fasta -r bam_alignment_qc.pdf ../../../01_processeddata/dummy_pool/align_calibration_strands/lambda_alignment_sorted.bam"
#     "-x	Do not plot per-reference information. Default: False"
# "bam_count_reads.py -z ../../../00_rawdata/reference_sequences/lambda/lambda_3.6kb.fasta -g -t bam_count_reads.tsv ../../../01_processeddata/dummy_pool/align_calibration_strands/lambda_alignment_sorted.bam"
# "bam_gc_vs_qual.py -f ../../../00_rawdata/reference_sequences/lambda/lambda_3.6kb.fasta -r bam_gc_vs_qual.pdf -t bam_gc_vs_qual.tsv ../../../01_processeddata/dummy_pool/align_calibration_strands/lambda_alignment_sorted.bam"
# "bam_multi_qc -h"
# "bam_ref_base_coverage.py -f ../../../00_rawdata/reference_sequences/lambda/lambda_3.6kb.fasta -t bam_ref_base_coverage.tsv ../../../01_processeddata/dummy_pool/align_calibration_strands/lambda_alignment_sorted.bam"
# "bam_soft_clips_tab.py -t bam_soft_clips_tab.tsv ../../../01_processeddata/dummy_pool/align_calibration_strands/lambda_alignment_sorted.bam"
# "bias_explorer.py -r bias_explorer.pdf bam_count_reads.tsv"

# EPI2ME:
# ~/_Temp/epi2me-cli-linux-2019.11.11-2920621 --dryrun -b 71622089 --consentedHuman -w 1760 --params '{"1":{"1_1_min_qscore":"7","1_2_detect_barcode":"RAB204"}}' \
# --inputfolder="/ag-halama/Microbiome/16S_Metabarcoding/01_processeddata/16S_nanopore_primer/Testpool/20191007_1559_MN31344_FAK76605_2bf006ff/basecall/pass/barcode01" \
# --outputfolder="/ag-halama/Microbiome/16S_Metabarcoding/02_analysis/16S_nanopore_primer/Testpool/20191007_1559_MN31344_FAK76605_2bf006ff/EPI2ME/barcode01"
# 
# When used with the --instance=ARG and -w options, it will download the summary of the workflow results as used in the report), e.g.
# --instance 12345 -t -w 1711
