## SETUP ##

## dependencies
import pandas as pd
import subprocess
import re
import matplotlib

## logging
sys.stdout = open(snakemake.log[0], 'w')
sys.stderr = open(snakemake.log[0], 'w')

## input files
input_dict = {
    'sam' : snakemake.input['sam'],
    'kronataxlist' : snakemake.input['kronataxlist'],
    'kronaseq2tax' : snakemake.input['kronaseq2tax']
}

## output files
output_dict = {
    'counttaxlist' : snakemake.output['counttaxlist']
}

# input_dict = {
#     'sam' : "01_processed_data/03_alignment/20191007_1559_MN31344_FAK76605_2bf006ff/barcode03/silva/filtered.sam",
#     'kronataxlist' : "METADATA/Reference_Sequences/silva/krona/species/taxlist.txt",
#     'kronaseq2tax' : "METADATA/Reference_Sequences/silva/krona/species/seqid2taxid.map"
# }
# output_dict = {
#     'counttaxlist' : "01_processed_data/03_alignment/20191007_1559_MN31344_FAK76605_2bf006ff/barcode03/silva_species/filtered.counttaxlist"
# }


## LOAD DATA ##

## list with path and taxID
df_taxlist = pd.read_csv(input_dict['kronataxlist'], sep='\t',
names=['pathname_slv','taxID_slv','rank_slv'],
usecols=['pathname_slv','taxID_slv'])

## list with accID and taxID
df_seq2tax = pd.read_csv(input_dict['kronaseq2tax'], sep='\t',
names=['accID','taxID_slv'])

## number of sam header lines
n_sam_header = int(subprocess.check_output(("awk 'BEGIN{count=0}; $1~/^@/{count++}; $1~!/^@/{exit}; END{print count}' " + input_dict['sam']),shell=True).decode("utf-8").rstrip())
## results from alignment
df_sam = pd.read_csv(input_dict['sam'], sep='\t',
header=None, usecols=[0, 2, 3, 4, 5, 9, 11, 12, 19],
#names=['QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL','NM','ms','AS','nn','tp','cm','s1','s2','de','SA/rl','rl/empty'],
#usecols=['QNAME','RNAME','POS','MAPQ','CIGAR','SEQ','NM','ms','de'], #!# sam files have differnent numbers of columns depending on SA, rl and zd flag
skiprows=n_sam_header)
df_sam.columns=['QNAME','RNAME','POS','MAPQ','CIGAR','SEQ','NM','ms','de']


## FILTERING READS ##

## get length of aligned read segment w/o softclips
cigar_len = []
for cigar in df_sam['CIGAR']:
    split_cigar = list(filter(None,re.split( "([0-9]+[A-Z])", cigar )))
    cigar_part = []
    for part in split_cigar:
        if re.search("[MI=X]",part):
            cigar_part.append(int(part[:-1]))
    cigar_len.append(sum(cigar_part))

## add to dataframe
df_sam['cigar_len'] = cigar_len

## convert "Gap-compressed per-base sequence divergence" into number
df_sam['de'] = df_sam['de'].str.extract(r'de:f:([\d.]+)', expand=False).astype('float')
#seems broken# df_sam['de'] = [float(mapq[6:]) for mapq in df_sam['de']]

## filter by aligned segment length
df_sam = df_sam.loc[((df_sam['cigar_len'] > int(snakemake.config['filtering']['len_min'])) &
(df_sam['cigar_len'] < int(snakemake.config['filtering']['len_max']))),:]
# df_sam = df_sam.loc[((df_sam['cigar_len'] > 1000) & (df_sam['cigar_len'] < 2800)),:]

## filter by Gap-compressed per-base sequence divergence
df_sam = df_sam.loc[(df_sam['de'] < 0.1),:]


## COUNTING READS ##

## get counts per accID (and convert to data frame)
counts_acc = pd.DataFrame(df_sam['RNAME'].value_counts()).reset_index().rename(columns={"index": "accID", "RNAME": "reads_here"})

## merge with meta data
counts_acc = pd.merge(counts_acc, df_seq2tax, how='left', on='accID').sort_values('taxID_slv')

## use multiindex to collapse counts per taxID
counts_acc.set_index(['taxID_slv', 'accID'], inplace=True)
df_combi = counts_acc.sum(level='taxID_slv').reset_index()

## filter reads below minimum required coverage
df_combi = df_combi.loc[(df_combi['reads_here'] >= int(snakemake.config['filtering']['min_featurereads'])),:]
# df_combi = df_combi.loc[(df_combi['reads_here'] >= 3),:]

## add path information
df_combi = pd.merge(df_combi, df_taxlist, how='left', on='taxID_slv')


## WRITE DATA ##

## open file for writing
wo = open(output_dict['counttaxlist'], "w")

## loop over dataframe rows
for idx, row in df_combi.loc[:,['reads_here','pathname_slv']].iterrows():
    ## split path into taxa
    pathname_split = row['pathname_slv'].split(';')
    null = pathname_split.pop(-1)
    ## combine to string with counts per taxon
    pathname_split = str(row['reads_here']) + "\t" + "\t".join(pathname_split) + "\n"
    null = wo.write(pathname_split)

## close file
wo.close()
