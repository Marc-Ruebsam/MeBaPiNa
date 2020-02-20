## SETUP ##

## dependencies
import pandas as pd

## logging
sys.stdout = open(snakemake.log[0], 'w')
sys.stderr = open(snakemake.log[0], 'w')

## input files
input_dict = {
    'report' : snakemake.input['report'],
    'taxlist' : snakemake.input['taxlist'] + "/taxlist.txt"
}

## output files
output_dict = {
    'ktaxlist' : snakemake.output['ktaxlist']
}

# input_dict = {
#     'report' : "01_processed_data/03_kmer_mapping/20191007_1559_MN31344_FAK76605_2bf006ff/barcode03/silva_species/filtered.kreport2",
#     'taxlist' : "METADATA/Reference_Sequences/silva/krona/species/taxlist.txt"
# }
# output_dict = {
#     'ktaxlist' : "01_processed_data/03_kmer_mapping/20191007_1559_MN31344_FAK76605_2bf006ff/barcode03/silva_species/filtered.ktaxlist"
# }

## LOAD DATA ##

## list with path and taxID
df_taxlist = pd.read_csv(input_dict['taxlist'], sep='\t', 
names=['pathname_slv','taxID_slv','rank_slv'], 
usecols=['pathname_slv','taxID_slv'])

## results from taxonomic classification
df_kreport = pd.read_csv(input_dict['report'], sep='\t', 
names=['percent','reads_cum','reads_here','rank_krak','taxID_slv','name'], 
usecols=['reads_here','taxID_slv'])
## exclude unclassified reads
df_kreport = df_kreport.loc[df_kreport['taxID_slv'] != 0,:]

## PROCESS DATA ##

## merge data frames
df_combi = pd.merge(df_kreport, df_taxlist, how='left', on='taxID_slv')

## open file for writing
wo = open(output_dict['ktaxlist'], "w")

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
