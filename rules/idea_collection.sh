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
# awk 'BEGIN{prnt_flag=1}; /fastq.gz/{next}; /Chromosome stats/{prnt_flag=0}; /summary section/{prnt_flag=1}; prnt_flag==1' qualimapReport.html > test.html

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
