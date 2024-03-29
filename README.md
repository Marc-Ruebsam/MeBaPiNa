created     27.01.2020
last edit:  21.06.2020 by Marc Ruebsam

================================
 Information text about the folder, what data is stored here and who is directly responsible for the contents and data structure.
================================

Project Folder: 16S_Metabarcoding
Description:    Microbiota profiling with long amplicons using nanopore sequencing:
                This folder is used for 16S Metabarcoding analysis using the Nanopore MinION platform.
                It contains subfolders for storage of raw data (created during the sequencing run), intermediate processing steps,
                analysis results and reports (created by the pipeline).
                The first three are organized per run (unique run ID), the last one is separated by time points and sample IDs.
                Additional information regarding the origin and processing of each sample can be received from the METADATA directory.
                The pipeline used for analysis of the data is based on the snakemake files in the Pipeline directory.
                Documentation contains written evaluation of the results, presentations and references.
Type of Data:   Fast5/Fastq/Pictures/Report files
Participants:   Silke Grauling-Halama (silke.grauling-halama@nct-heidelberg.de), Pornpimol Charoentong (Pornpimol.Charoentong@nct-heidelberg.de), Marc Ruebsam (marc.ruebsam@bioquant.uni-heidelberg.de)
Correspondence: Silke Grauling-Halama (silke.grauling-halama@nct-heidelberg.de)


================================
 Data creation and handling
================================

The raw data is created by the MinKNOW software during the MinION run.
Crucial information regarding the samples, library preparation and sequencing run are stored inside METADATA/EXPERIMENT_SEQUENCING.xlsx file.

== Sequencing with MinKNOW ==

- Follow the protocol for flow cell preparation and loading of the library
    - Note down the number of pores during the flow cell check
- Select the folder specified as Project (see above) as output directory (will always be a "00_raw_data" folder inside the project specific folder)
- Specify a experiment and sample name (as above)
    - Please do not use spaces or special characters (e.g. äöüß.,<>|), except "_"
- Activate VBZ compression
- Deactivate live basecalling
- Start the sequencing run
- MinKNOW will create the output directory inside 00_raw_data according to the experiment and sample names specified
    - if you want to create and save PDF reports, create a folder inside the Run ID folder called "reports" (the same folder that contains the fast5 folder)
    - move PDFs in there

== After the run has finished ==

- Make sure the sequencing run has finished successfully and the raw data is in the expected location
- Close all open folders and files related to the current run

================================
 Analysis with MeBaPiNa Pipeline
================================

MeBaPiNa (metabarcoding analysis pipeline for Nanopore datasets) is a pipeline implemented in snakemake.
It takes raw fast5 files and automatically processes it according to the specifications in the METADATA/PIPELINE_CONFIG.yaml file.
Statistics and figures for the requested samples are reported in the 03_report directory.

== Installing snakemake ==

- Has to be done once on the machine used for running the Pipeline
    - PLEASE THINK BEFORE COPY & PASTING
    - Install conda
        ## update machine
        sudo yum update -y
        sudo yum upgrade -y
        sudo yum install -y epel-release ## EPEL repository for yum
        sudo yum install -y cifs-utils ## mounting of windows network shares
        sudo yum install -y wget ## querry downloads
        sudo yum install -y git ## git
        ## get lateste miniconda
        wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
        ## create directory for miniconda owed by user
        sudo mkdir -p /opt/miniconda
        pid_me=$(id --user)
        gid_me=$(id --group)
        sudo chown ${pid_me}:${gid_me} /opt/miniconda
        unset pid_me gid_me
        ## execute installation script
        bash miniconda.sh -b -p /opt/miniconda/miniconda3
        export PATH="/opt/miniconda/miniconda3/bin:$PATH"
        rm miniconda.sh
        ## update and initialize
        hash -r
        conda update -y -q conda
        conda init bash
    - Install snakemake (create snakemake environment with conda)
        ## create conda environment of snakemake
        conda create -y -c defaults -c bioconda -c conda-forge -n snakemake snakemake xlrd
                ## alternatively use the configuration for 5.10
                conda env create -f Pipeline/MeBaPiNa/envs/snakemake.yml -n snakemake
        ## remove temp files
        conda clean -y --all
        - Install guppy standalone
                - download guppy from
                https://europe.oxfordnanoportal.com/software/analysis/ont_guppy_3.4.1-1~xenial_amd64.deb
                - install

== Before running the Pipeline ==

- Update the METADATA/PIPELINE_CONFIG.yaml file
    - Experiment specifications
        - project, tmp:             build the path to the project directory like {project}/{tmp}
        - meta:                     path to the EXPERIMENT_SEQUENCING.xlsx inside the project directory like {project}/{tmp}/{meta}
        - log:                      path to the ANALYSIS_PROGRESS_MANAGEMENT.csv inside the project directory like {project}/{tmp}/{log}
        - samples:                  list of sample names from the EXPERIMENT_SEQUENCING.xlsx to analze
    - Methods for analysis
        - methodologie:             list of method abbreviations to use for analysis
    - Workstation specifications
        - gpu:                      whether or not a guppy compatible GPU is available
    - Filtering options
        - q_min:                    minimal quality score for read filtering
        - len_min:                  minimal length for read filtering
        - len_max:                  maximal length for read filtering
        - min_featurereads          minimum number of reads per cluster or taxon
        - min_readidentity          minimum read identity for clustering (OTU)
        - min_confidence            confidence used for classification of taxa at {refrank} level
        - plot_sample               down-sampling reads in plots to this number to reduce computational complexity (distorts statistics shown in some plots)
    - Reference database (currently only {refsource} is implemented)
        - source:                   which reference database to use
        - rank:                     at which rank the taxonomy should be determined

== Run MeBaPiNa ==

I. Run full pipeline

- Activate environment if used conda is used for snakemake
        ## activate conda environment
        conda activate snakemake
- Change directory into the project directory (or use absolute paths below)
- DRY RUN: Optionally execute the pipeline as a dry run showing all rules that are going to be executed
    - the '-n' invokes a dry run and '-pr' includes higher verbosity
    - change the path to the config file and Pipeline snakefile if necessary
        ## dry run pipeline with verbosity
        snakemake -npr --use-conda --cores 'all' --configfile METADATA/PIPELINE_CONFIG.yaml --snakefile Pipeline/MeBaPiNa/Snakefile
- Execute the pipeline
    - change the path to the config file and Pipeline snakefile if necessary
        ## run pipeline
        snakemake -pr --use-conda --cores 'all' --configfile METADATA/PIPELINE_CONFIG.yaml --snakefile Pipeline/MeBaPiNa/Snakefile
- Please note that the --cores 'all' might not work to request all available threads at the machine and has to be changed to an integer manually

II. Perform only basecalling

- In case you only want to basecall the data and have a look at the data before running the full pipeline
- Or in case you want to perform basecalling on a separate machine than the rest of the pipeline
    - you can run
        ## run the pipeline only until basecalling has been finished
        snakemake -pr --use-conda --cores 'all' --configfile METADATA/PIPELINE_CONFIG.yaml --snakefile Pipeline/MeBaPiNa/Snakefile only_basecall
    - output will be generated inside the 03_report/ dorectory and ccan be reviewed
    - please note that this instance of the pipeline has to finish before you can start further analysis by invoking the command from I.
    - in other words: no dependecy between this instance and another can be established automatically and has to be insured manually

III. Update report file / add missing reports

- In case you want to update the log file ANALYSIS_PROGRESS_MANAGEMENT.csv
    - you can run
        ## rerun report rules and any of its upwards dependencies
        snakemake -pr --use-conda --cores 'all' --configfile METADATA/PIPELINE_CONFIG.yaml --snakefile Pipeline/MeBaPiNa/Snakefile update_report

IV. Rerun certain step

- You can rerun a certain step in the pipeline or reproduce a certain file (e.g. if it was corrupted)
  - find the rule corresponding to the processing step in question
    - check the output location for the output files
    - if the output files still exist, delete or move them (rerunning only works if the file doesn't already exist)
  - get the path to the desired output file and call snakemake with it, e.g.:
      ## rerun krona plot for one of the samples
      snakemake -pr --use-conda --cores 'all' --configfile METADATA/PIPELINE_CONFIG.yaml --snakefile Pipeline/MeBaPiNa/Snakefile 16S_Metabarcoding/02_analysis_results/03_kmer_mapping/20191007_1559_MN31344_FAK76605_2bf006ff/barcode06/silva_Species/krona.html
  - you can list multiple files to rerun all of them
- Note! rerunning of an intermediate step does not automatically update all below dependencies
  - if this is desired, proceed as above, by deleting or moving the file in question
  - then call the pipeline as normal

== After running the pipeline ==

- Check the snakemake log file (in .snakemake/logs/) for any errors
- Repeat a dry run to check that only the "all_target" rule is left (this rule will be executed every time)
- check the results inside the 03_report directory


================================
 Directory contents and pipeline processing steps
================================

The following list gives an overview of the files inside the 16S_Metabarcoding directory.
The list does not include temporary files or log files.
Please have a look at the Documentation/ThesisManuscript/Figures/DAG.pdf for an overview of the rule dependencies.

== METADATA ==

./METADATA/
  - directory containing important metadata files and reference sequences

    PIPELINE_CONFIG.yaml
      - configurations of the pipeline
      - see above

    Reference_Sequences/
      - directory containing reference sequences

        primers/
          - primers used for the extraction of the amplicon region from the references

        {refsource}/
          - directory containing the reference database files and subdirectories for the methodology specific files

            ## rule download_reffiles and construct_refseq:
            reference.fasta
              - reference sequences trimmed to amplicon region, filtered by length and removed duplicates (100% identity)

            ## rule indexing_reference:
            reference.mmi
              - reference file for alignments with minimap2

                ## rules construct_reftax and building_database:
                kraken2/
                  - kraken2 and bracken specific reference database files

                ## rules construct_reftax, q2import_reftax, q2import_refseq and q2train_classifyer:
                qiime/
                  - qiime specific reference database files
                  - these files might be incomplete due to an unsuccessful construction

                ## rules construct_reftax:
                krona/
                  - krona specific reference database files

        lambda/
          - lambda phage reference sequence

        zymobiomics/
          - reference files only containing only species present in the zymobiomics mock community

== Pipeline ==

./Pipeline/
  - a directory containing software and scripts associated with the pipeline

    MeBaPiNa/
      - directory containing files required to run the pipeline

        Snakefile
          - central file of the pipeline
          - starts execution and sources other files

        envs/
          - conda environment files used in the pipeline

        rules/
          - definition of rules used in the pipeline

        scripts/
          - scripts used in the rules
          - some additional scripts

== Raw data ==

./00_raw_data/
  - directory containing raw sequencing files of all runs
  - priority:
    - these files are of highest priority and cannot be replaced

    ## rule move_raw:
    {run}/
      - run specific sub-directory

        fast5/
          - directory containing raw fast5 sequencing files

        reports/
          - optional directory with run reports as PDFs
          - created manually in MinKNOW (feature automation was promised by ONT)

== Intermediate processing files ==

./01_processed_data/
  - directory containing intermediate processing files of all steps and runs
  - priority:
    - these files are produced during key steps in the analysis processes
    - they can be replaced if the files from ./00_raw_data/ are still present, but this might require a considerable amount of recomputation
    - they are required in case some steps should be reran (e.g. with new parameters)

    01_basecalling/
      - directory containing the basecalled raw-reads of all runs

        {run}/
          - run specific sub-directory

            ## rule basecall_raw:
            guppy_basecaller_logs/
              - logs from basecaller
            sequencing_summary.txt
              - table with information about all basecalled reads
            pass/
              - directory containing basecalled fastq files of all barcodes
                {barcode}
                  - sub-directories with fastq files assigned to the barcode
                unclassified
                  - sub-directory with fastq files without detected barcode (contains all files if library wasn't multiplexed)

            ## rule sorting_seqsum_barc and splitting_seqsum_barc
            sequencing_summary/
              - same information as in sequencing_summary.txt
                sequencing_summary_sorted.txt
                  - but sorted by barcode column and start time
                split/
                  - and as separate files per barcode (can be used to extract raw reads per barcode)

    02_trimming_filtering/
      - directory containing the trimmed and filtered basecalled-reads of all runs

        {run}/
          - run specific sub-directory

            {barcode}/
              - barcode specific sub-directory

                ## rule trim_basecalled
                trimmed.fastq
                  - trimming and second demultiplexing of basecalled reads
                other/
                  - reads assigned to other barcodes during second demultiplexing
                  - these reads are not used for anything and could be deleted

                ## rule filter_trimmed:
                filtered.fastq
                  - length and quality filtered filtered reads

    03_alignment/
      - directory containing the alignment of trimmed and filtered reads of all runs

        {run}/
          - run specific sub-directory

            {barcode}/
              - barcode specific sub-directory

                {refsource}/
                  - reference database specific sub-directory

                    ## rule aligning_filtered and filter_aligned:
                    filteredsorted.bam
                      - filtered alignement of the trimmed and filtered reads
                    filteredsorted.bam.bai
                      - index of bam

    03_kmer_mapping/
      - directory containing the k-mer mapping of trimmed and filtered reads of all runs

        {run}/
          - run specific sub-directory

            {barcode}/
              - barcode specific sub-directory

                {refsource}_{refrank}/
                  - reference database and rank specific sub-directory

                    ## rule kmermap_filtered:
                    filtered.kraken2
                      - per read classification information
                    filtered.kreport2
                      - per taxa classification information

                    ## rule retax_kmermap:
                    Species.bracken
                      - per taxa reestimation information
                    Species.kreport2
                      - per taxa classification information after reestimation

    03_otu_picking/
      - directory containing the clustered OTUs and the taxonomic assignmnts of trimmed and filtered reads of all runs

        {run}/
          - run specific sub-directory

            {barcode}/
              - barcode specific sub-directory

                {refsource}/
                  - reference database specific sub-directory

                    ## rule generate_manifestfile, q2import_filtered, q2derep_imported, q2otupick
                    cluster_centseq.qza
                      - center sequence per OTU of the trimmed and filtered reads
                    cluster_ftable.qza
                      - table with counts per OTU of the trimmed and filtered reads
                    cluster_newrefseq.qza
                      - new reference sequences
                      - can be used for subsequent runs of open-reference clustering for consistent definitions of features across open-reference feature tables.

                    ## q2uchime_otus and q2filter_uchime:
                    filt_ftable.qza
                      - table with counts per OTU after chimera removal and min coverage filtering
                    filt_centseq.qza
                      - center sequence per OTU after chimera removal and min coverage filtering

                    {other directories}
                      - storage location for temporary files during conversion

                {refsource}_{refrank}/
                  - reference database and rank specific sub-directory

                    ## rule rereplicate_q2filter and kmermap_q2rereplicate:
                    filtered.kraken2
                      - per OTU classification information
                    filtered.kreport2
                      - per taxa classification information

== Result files ==

./02_analysis_results/
  - directory containing result files of all steps and runs
  - priority:
    - these files are of low priority
    - they can easily be replaced if the files from ./01_processed_data/ are still present
    - important files are copied to ./03_report/

    01_basecalling/
      - directory containing plots and statistics of the basecalled raw-reads of all runs

        {run}/
          - run specific sub-directory

            ## rule plot_fastq_pipe_basecall:
            fastqc/
              - read QC: all passed reads

            ## rule plot_nanocomp_seqsum_basecall:
            nanocomp/
              - barcode QC: per barcode

            ## rule plot_nanoplot_seqsum_basecall:
            nanoplot/
              - general QC: all reads, including calibtation strads

            ## rule plot_nanoqc_fastq_basecall:
            nanoqc/
              - per base QC: all reads

            ## rule plot_pycoqc_seqsum_basecall:
            pycoqc/
              - general QC: all reads

    02_trimming_filtering/
      - directory containing plots and statistics of the trimmed and filtered basecalled-reads of all runs

        {run}/
          - run specific sub-directory

            ## rule plot_fastqc_fastq_filter:
            fastqc/
              - read QC: trimmed and filtered barcoded reads

            ## rule plot_nanocomp_fastq_filter:
            nanocomp/
              - barcode QC: trimmed and filtered barcoded reads

            ## rule plot_nanoplot_fastq_filter:
            nanoplot/
              - general QC: trimmed and filtered barcoded reads

            ## rule plot_fastq_pipe_filter:
            nanoqc/
              - per base QC: trimmed and filtered barcoded reads

    03_alignment/
      - directory containing plots and statistics of the alignment of trimmed and filtered reads of all runs

        {run}/
          - run specific sub-directory

            {barcode}/
              - barcode specific sub-directory

                {refsource}/
                  - reference database specific sub-directory

                      ## rule plot_refseq_coverage:
                      covdist.pdf
                      covpos.pdf
                        - coverage of reference sequences
                        - covdist.pdf seem broken when executing on the VM

                      ## rule plot_pycoqc_aligned:
                      pycoqc.html
                      pycoqc.json
                        - alignment QC

                {refsource}_{refrank}/
                  - reference database and rank specific sub-directory

                      ## rule counttax_aligned:
                      aligned.counttaxlist
                        - taxonomic classification

                      ## rule plot_krona_aligned_text:
                      krona.html
                        - visualization of taxonomic classification

    03_kmer_mapping/
      - directory containing plots and statistics of the k-mer mapping of trimmed and filtered reads of all runs

        {run}/
          - run specific sub-directory

            {barcode}/
              - barcode specific sub-directory

                {refsource}_{refrank}/
                  - reference database and rank specific sub-directory

                      ## rule counttax_kmermap:
                      kmer.counttaxlist
                        - taxonomic classification

                      ## rule plot_krona_kmermap_kraken:
                      krona.html
                        - visualization of taxonomic classification

                      ## rule plot_krona_kmermap_text:
                      krona_bracken.html
                        - visualization of taxonomic classification after reestimation

    03_otu_picking/
      - directory containing plots and statistics of the clustered OTUs and the taxonomic assignmnts of trimmed and filtered reads of all runs

        {run}/
          - run specific sub-directory

            {barcode}/
              - barcode specific sub-directory

                {refsource}/
                  - reference database specific sub-directory

                      ## rule plot_qiime2_q2otupick:
                      q2otupick\
                        - clustered reads overview

                      ## rule plot_qiime2_q2filter:
                      q2filter\
                        - clustered reads overview after filtering

                {refsource}_{refrank}/
                  - reference database and rank specific sub-directory

                      ## rule counttax_q2kmermap:
                      kmer.counttaxlist
                        - taxonomic classification

                      ## rule plot_krona_q2rerep:
                      krona.html
                        - visualization of taxonomic classification

== Report per sample ==

./03_report/
  - directory containing result files, plots and statistics per timepoint and sample
  - priority:
    - the files in this directory are of high priority
    - they can be recreated if the files from ./01_processed_data/ and ./02_analysis_results/ are still presentat
    - they represent the final output of the pipeline

    Reference_Sequences/
      - directory with plots and statistics of reference database

        {refsource}/
          - reference database specific sub-directory

            ## rule stat_refseq_lenstat:
            reference_lengthdist.pdf
            reference_lengthdist.tsv
              - length of reference sequences

            ## rule stat_refseq_taxaranks
            reference_taxaranks.tsv
              - distribution of reference taxa ranks

    {timepoint}/
      - timepoint specific folder or non-PROMISE

        {sample}/
          - sample specific folder

            {run}-{barcode}/
              - sample per run specific folder
              - in case the same sample is sequenced in multiple runs

                ## rule stat_general_readbasecount:
                read_base_counts.tsv
                  - statistics of reads and bases of raw, basecalled, trimmed and filtered reads

                ## rule :
                03_alignment-{refsource}-alignment_rates.tsv
                  - aligment statistics and error rates

                ## rule stat_otu_feature:
                03_otu_picking-{refsource}-feature_counts.tsv
                  - number of clusters and reads before and after filtering

                ## rule stat_align_taxa:
                03_alignment-{refsource}_{refrank}-taxa_counts.tsv
                  - statistics of classified taxa and reads

                ## rule stat_kmer_taxa:
                03_kmer_mapping-{refsource}_{refrank}-taxa_counts.tsv
                - statistics of classified taxa and reads

                ## rule stat_kmer_retaxa:
                03_kmer_mapping-{refsource}_{refrank}-retaxa_counts.tsv
                - statistics of classified taxa and reads after reestimation

                ## rule stat_otu_taxa:
                03_otu_picking-{refsource}_{refrank}-taxa_counts.tsv
                - statistics of classified taxa and reads

                ## rule stat_align_taxa_diversity:
                03_alignment-{refsource}_{refrank}-taxa_covdist.pdf
                - distribution of taxa abundance
                03_alignment-{refsource}_{refrank}-taxa_diversity.tsv
                - diversity, richness and evenness of community

                ## rule stat_kmer_taxa_diversity:
                03_kmer_mapping-{refsource}_{refrank}-taxa_covdist.pdf
                  - distribution of taxa abundance
                03_kmer_mapping-{refsource}_{refrank}-taxa_diversity.tsv
                  - diversity, richness and evenness of community

                ## rule stat_otu_taxa_diversity:
                03_otu_picking-{refsource}_{refrank}-taxa_diversity.tsv
                - distribution of taxa abundance
                03_otu_picking-{refsource}_{refrank}-taxa_covdist.pdf
                - diversity, richness and evenness of community

                {other files}
                  - are copies from the similarly named files in the ./02_analysis_results/ directory
