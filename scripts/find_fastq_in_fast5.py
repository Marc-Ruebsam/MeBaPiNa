import os
import sys
import pandas as pd
import subprocess
from datetime import datetime
from Bio import SeqIO
from statistics import mean
from math import log10
from gzip import open as gopen
from re import split as rsplit

sys.stdout = open(snakemake.log[0], 'w')
sys.stderr = open(snakemake.log[0], 'w')

###########################
## EXPERIMENT START TIME ##
###########################

## get experiment start from sequencing_telemetry.js
exp_start_time = str(subprocess.check_output("grep exp_start_time " + 
os.path.join(*snakemake.input[0].split("/")[0:-1], "sequencing_telemetry.js") + 
" | head -n 1", shell = True)).split('"')[3]

## convert to datetime
exp_start_time = datetime.strptime(exp_start_time, "%Y-%m-%dT%XZ")

################
## FASTQ DATA ##
################

## initialize dataframe
with open(snakemake.output[0], 'w') as wo:
    wo.write( "\t".join((
    "filename_fastq",
    "read_id",
    "run_id",
    "channel",
    "start_time",
    "sequence_length_template",
    "mean_qscore_template",
    "barcode_arrangement" )) + "\n" )

    ## get list of all files (and directories) in pass directory
    for root, dirs, files in os.walk(snakemake.input[0], topdown=False):
        
        ## loop over all files
        for fastq_file in files:
            ## get file from list
            fastq_file_path = os.path.join(root, fastq_file)
            
            ## open file
            with gopen(fastq_file_path, 'rt') as f:
                
                ## loop over sequences
                for l in SeqIO.parse(f, "fastq"):
                    ## file name
                    filename_fastq = fastq_file
                    ## barcode
                    barcode_arrangement = root.split("/")[-1]
                    
                    ## information from read description
                    read_id, run_id, sampleid, read, channel, start_time = rsplit(" \w+=", l.description)
                    ## time since experiment start
                    start_time = datetime.strptime(start_time, "%Y-%m-%dT%XZ")
                    start_time = str((start_time - exp_start_time).total_seconds() + 3600)
                    ## qscore 
                    mean_qscore_template = mean([ 10**(i*(-0.1)) for i in l.letter_annotations["phred_quality"] ]) # round(mean(l.letter_annotations["phred_quality"]),6)
                    mean_qscore_template = str(round((-10)*log10(mean_qscore_template),6))
                    ## length
                    sequence_length_template = str(len(l.seq))
                    
                    ## create row in data
                    wo.write( "\t".join((
                    filename_fastq,
                    read_id,
                    run_id,
                    channel,
                    start_time,
                    sequence_length_template,
                    mean_qscore_template,
                    barcode_arrangement
                    )) + "\n" )
