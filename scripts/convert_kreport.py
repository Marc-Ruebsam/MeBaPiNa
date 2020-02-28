## SETUP ##

## dependencies
import pandas as pd

## logging
sys.stdout = open(snakemake.log[0], 'w')
sys.stderr = open(snakemake.log[0], 'w')

## input files
input_dict = {
    'kreport' : snakemake.input['kreport'],
    'kronataxlist' : snakemake.input['kronataxlist']
}

## output files
output_dict = {
    'counttaxlist' : snakemake.output['counttaxlist']
}

# input_dict = {
#     'kreport' : "01_processed_data/03_kmer_mapping/20191007_1559_MN31344_FAK76605_2bf006ff/barcode03/silva_species/filtered.kreport2", ## kraken2
#     # 'kreport' : "02_analysis_results/03_kmer_mapping/20191007_1559_MN31344_FAK76605_2bf006ff/barcode03/silva_species/Species.kreport2", ## bracken
#     'kronataxlist' : "METADATA/Reference_Sequences/silva/krona/species/taxlist.txt"
# }
# output_dict = {
#     'counttaxlist' : "02_analysis_results/03_kmer_mapping/20191007_1559_MN31344_FAK76605_2bf006ff/barcode03/silva_species/kmer.counttaxlist"
# }

## LOAD DATA ##

## list with path and taxID
df_taxlist = pd.read_csv(input_dict['kronataxlist'], sep='\t', 
names=['pathname_slv','taxID_slv','rank_slv'], 
usecols=['pathname_slv','taxID_slv'])

## results from taxonomic classification
df_kreport = pd.read_csv(input_dict['kreport'], sep='\t', 
names=['percent','reads_cum','reads_here','rank_krak','taxID_slv','name'], 
usecols=['reads_here','taxID_slv'])
## exclude unclassified reads
df_kreport = df_kreport.loc[df_kreport['reads_here'] != 0,:]

## PROCESS DATA ##

## merge data frames
df_combi = pd.merge(df_kreport, df_taxlist, how='left', on='taxID_slv')
## add "path" string to unclssified
df_combi.loc[df_combi['taxID_slv'] == 0,'pathname_slv'] = "unclassified;"

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
