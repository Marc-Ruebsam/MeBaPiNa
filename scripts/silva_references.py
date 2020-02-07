## SETUP ##

## dependencies
import pandas as pd

## input files
input_dict = {
# 'fasta' : "../SILVA_132_SSURef_Nr99_tax_silva.fasta",
# 'accID2taxid' : "SILVA_132_SSURef_Nr99_tax_silva/taxonomy/tax_slv_ssu_132.acc_taxid" ,
'taxlist' : "SILVA_132_SSURef_Nr99_tax_silva/taxonomy/tax_slv_ssu_132.txt",
'slvmap' : "SILVA_132_SSURef_Nr99_tax_silva/taxonomy/taxmap_slv_ssu_ref_nr_132.txt",
'ncbimap' : "SILVA_132_SSURef_Nr99_tax_silva/taxonomy/taxmap_embl_ssu_ref_nr99_132.txt",
'kraknodes' : "kraken2_silva_132/taxonomy/nodes.dmp"
}
# input_dict = {
# 'fasta' : "../SILVA_138_SSURef_NR99_tax_silva.fasta",
# 'taxlist' : "tax_slv_ssu_138.txt",
# 'accID2taxid' : "tax_slv_ssu_138.acc_taxid" ,
# 'ncbimap' : "taxmap_embl-ebi_ena_ssu_ref_138.txt"
# }
## output files
output_dict = {
'krona' : "krona_silva_132/krona.tab"
}

## IMPORTING TABLES ##

## load taxID list
df_taxlist = pd.read_csv(input_dict['taxlist'], sep='\t', 
names=['pathname_slv','taxID_slv','rank_slv','remark','release'], 
usecols=['pathname_slv','taxID_slv','rank_slv'])
df_taxlist['depth'] = [len(splt)-1 for splt in df_taxlist['pathname_slv'].str.split(';')]
df_taxlist['name'] = [splt[-2] for splt in df_taxlist['pathname_slv'].str.split(';')]

## load SILVA taxIDs w/o species classification
df_slvmap = pd.read_csv(input_dict['slvmap'], sep='\t', 
skiprows=1,
names=['accID','start','end','path_slv','name','taxID_slv'],
usecols=['accID','name','taxID_slv'])

## load NCBI taxIDs w/ species classification
df_ncbimap = pd.read_csv(input_dict['ncbimap'], sep='\t', 
# skiprows=1, ## needed for SILVA 138
names=['accID','start','end','path_ncbi','name','taxID_ncbi'],
usecols=['accID','start','end','path_ncbi','taxID_ncbi'])







## check if all used taxons (assigned to sequences) are in the txon list
if set(df_slvmap['taxID_slv']) - set(df_taxlist['taxID_slv']) != set():
    print("Warning: not all used taxons (" + input_dict['slvmap'] + ") are in the taxon list (" + input_dict['taxlist'] + ").")

## merge silva dataframe data frames
df_slv = pd.merge(df_slvmap, df_taxlist, how='left', on='taxID_slv')

## check if all used taxons (assigned to sequences) are in the taxID targetID assignment
if set(df_slv['taxID_slv']) - set(df_kraknodes['taxID_slv']) != set():
    print("Warning: not all used taxons (" + input_dict['slvmap'] + ") are in the taxID targetID assignment (" + df_kraknodes['taxlist'] + ").")

## merge silva dataframe data frames
df_slv = pd.merge(df_slv, df_kraknodes, how='left', on='taxID_slv')



if set(df_all['accID']) - set(df_accID2taxid['accID']) != set():
    print("Warning: accIDs form fasta file don't match taxmap_embl_ssu_ref_nr99_132.txt, excluding missing entries.")

df_all = pd.merge(df_all, df_ncbimap, how='inner', on='accID')


## CREATE NEW TAXONS ##

## create new species taxID by combining SILVA and ENA taxIDs (use different Ids, simple combination introduces overlap)
df_all['taxID_new'] = (df_all['taxID_ncbi'].map(str) + df_all['taxID_x'].map(str)).str.replace('-', '') ## some taxIDs in the ena file are negative oO
# if set(df_taxlist['taxID']).intersection(set(df_all['taxID_new'])) != set():
#     print("Warning: Have to double down on taxID concatenation because of overlap with SILVA taxIDs.")
#     df_all['taxID_new'] = (df_all['taxID_new'].map(str) + df_all['taxID_ncbi'].map(str)).map(int)

# enaIDs = wget "https://www.ebi.ac.uk/ena/data/view/AY929283&display=xml&download=xml&filename=AY929283.xml" -q -O -| awk -F "\"" '/taxId/{print $4}'
# enaIDs = wget "https://www.ebi.ac.uk/ena/data/view/Taxon:1578&display=xml&download=xml&filename=1578.xml" -q -O - | awk -F"=" '/^<taxon scientificName/{split($2,sp2,"\""); split($3,sp3,"\""); split($5,sp5,"\""); print sp2[2]"\t"sp3[2]"\t"sp5[2]}'

## list of possible ranks (I am not a microbiologist and this system is confusing)
tax_rank = pd.DataFrame( { "rank" : [ "superdomain",
"domain","major_clade","infradomain",
"superkingdom","kingdom","subkingdom","infrakingdom",
"superphylum","phylum","subphylum","infraphylum",
"superclass","class","subclass","infraclass",
"superorder","order","suborder","infraorder",
"superfamily","family","subfamily","infrafamily",
"supergenus","genus",
"name" ],
"index" : list(range(27)) } )
## assign new rank one below old rank
idx = pd.merge(df_all['rank'], tax_rank, how="left", on="rank")['index'] + 1
df_all['rank_new'] = tax_rank.loc[idx,'rank'].values


## WRITE FILES ##

## write new taxID list
df_taxlist.to_csv(
"tax_slv_ssu_132_species.txt", mode='w', sep='\t', header=False, index=False, quoting=0, float_format="%g")
pd.concat([
    df_all['path'].map(str) + df_all['name_y'].map(str) + ";",
    df_all['taxID_new'],
    df_all['rank_new'],
    pd.Series(["w"] * len(df_all.index)),
    pd.Series([0] * len(df_all.index))
], axis=1).drop_duplicates(subset=[0,1]).to_csv(
"tax_slv_ssu_132_species.txt", mode='a', sep='\t', header=False, index=False, quoting=0, float_format="%g")
## write new accID to taxID map
df_all.loc[:,['accID','taxID_new']].sort_values('accID').to_csv(
"tax_slv_ssu_132_species.acc_taxid", mode='w', sep='\t', header=False, index=False, quoting=0, float_format="%g")


































## load kraken2 taxID targetID assignment
df_kraknodes = pd.read_csv(input_dict['kraknodes'], sep='\t', 
names=['taxID_slv','line1','targetID_slv','line2','rank_slv','line3','remark','line4'], 
usecols=['taxID_slv','targetID_slv'])


## KRONA ##

## check if all kraken2 taxons are in the silva txon list
if set(df_taxlist['taxID_slv']) - set(df_kraknodes['taxID_slv']) != set():
    print("Warning: not all silva taxons (" + input_dict['taxlist'] + ") are represented by kraken2 taxon nodes (" + input_dict['kraknodes'] + ").")

## merge silva dataframe data frames
df_krona = pd.merge(df_taxlist, df_kraknodes, how='left', on='taxID_slv')

## select columns
df_krona = df_krona.loc[:,['taxID_slv','depth','targetID_slv','rank_slv','name']]

## remove useless names
df_krona.loc[df_krona['name'].isin(["uncultured","Incertae Sedis","Unknown Family"]),'rank_slv'] = "no rank"

## add root
df_krona = df_krona.append(pd.Series([1,0,1,"no rank","root"], index=['taxID_slv','depth','targetID_slv','rank_slv','name']), ignore_index=True)

## write new taxID list
df_krona.to_csv(
output_dict['krona'], mode='w', sep='\t', header=False, index=False, quoting=0, float_format="%g")
