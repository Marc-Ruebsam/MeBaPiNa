## SETUP ##

## dependencies
import pandas as pd
import os

## logging
sys.stdout = open(snakemake.log[0], 'w')
sys.stderr = open(snakemake.log[0], 'w')

## input files
input_dict = {
    'taxlist' : snakemake.input['taxlist'],
    'slvmap' : snakemake.input['slvmap'],
    'ncbimap' : snakemake.input['ncbimap'],
    'ncbikrona' : snakemake.input['ncbikrona']
}

## output files
output_dict = {
    'kraknames_S' : snakemake.output['kraknames_S'],
    'kraknodes_S' : snakemake.output['kraknodes_S'],
    'krakseq2tax_S' : snakemake.output['krakseq2tax_S'],
    'kraknames_G' : snakemake.output['kraknames_G'],
    'kraknodes_G' : snakemake.output['kraknodes_G'],
    'krakseq2tax_G' : snakemake.output['krakseq2tax_G'],
    'kronataxtab_S' : snakemake.output['kronataxtab_S'],
    'kronataxlist_S' : snakemake.output['kronataxlist_S'],
    'kronaseq2tax_S' : snakemake.output['kronaseq2tax_S'],
    'kronataxtab_G' : snakemake.output['kronataxtab_G'],
    'kronataxlist_G' : snakemake.output['kronataxlist_G'],
    'kronaseq2tax_G' : snakemake.output['kronaseq2tax_G']
}

# input_dict = {
#     'ncbikrona' : "METADATA/Reference_Sequences/ncbi/krona/taxonomy.tab",
#     'ncbimap' : "METADATA/Reference_Sequences/silva/ncbimap.txt",
#     'slvmap' : "METADATA/Reference_Sequences/silva/slvmap.txt",
#     'taxlist' : "METADATA/Reference_Sequences/silva/taxlist.txt"
# }
# output_dict = {
#     'kraknames_S' : "METADATA/Reference_Sequences/silva/kraken2/species/taxonomy/names.dmp",
#     'kraknodes_S' : "METADATA/Reference_Sequences/silva/kraken2/species/taxonomy/nodes.dmp",
#     'krakseq2tax_S' : "METADATA/Reference_Sequences/silva/kraken2/species/seqid2taxid.map",
#     'kraknames_G' : "METADATA/Reference_Sequences/silva/kraken2/genus/taxonomy/names.dmp",
#     'kraknodes_G' : "METADATA/Reference_Sequences/silva/kraken2/genus/taxonomy/nodes.dmp",
#     'krakseq2tax_G' : "METADATA/Reference_Sequences/silva/kraken2/genus/seqid2taxid.map",
#     'kronataxtab_S' : "METADATA/Reference_Sequences/silva/krona/species/taxonomy.tab",
#     'kronataxtab_G' : "METADATA/Reference_Sequences/silva/krona/genus/taxonomy.tab",
#     'kronataxlist_S' : "METADATA/Reference_Sequences/silva/krona/species/taxlist.txt",
#     'kronataxlist_G' : "METADATA/Reference_Sequences/silva/krona/genus/taxlist.txt",
#     'kronaseq2tax_S' : "METADATA/Reference_Sequences/silva/krona/species/seqid2taxid.map",
#     'kronaseq2tax_G' : "METADATA/Reference_Sequences/silva/krona/genus/seqid2taxid.map"
# }
# krak_dir_S = "/".join(output_dict['krakseq2tax_S'].split("/")[:-1])
# krak_dir_G = "/".join(output_dict['krakseq2tax_G'].split("/")[:-1])
# krona_dir_S = "/".join(output_dict['kronataxtab_S'].split("/")[:-1])
# krona_dir_G = "/".join(output_dict['kronataxtab_G'].split("/")[:-1])
# os.makedirs( krak_dir_S + "/library" )
# os.makedirs( krak_dir_S + "/taxonomy" )
# os.makedirs( krak_dir_G + "/library" )
# os.makedirs( krak_dir_G + "/taxonomy" )
# os.makedirs( krona_dir_S )
# os.makedirs( krona_dir_G )


## LOAD DATA ##

## load taxID list
df_taxlist = pd.read_csv(input_dict['taxlist'], sep='\t', 
names=['pathname_slv','taxID_slv','rank_slv','remark','release'], 
usecols=['pathname_slv','taxID_slv','rank_slv'])

## load SILVA taxIDs w/o species classification
df_slvmap = pd.read_csv(input_dict['slvmap'], sep='\t', 
skiprows=1,
names=['accID','start','end','path_slv','name','taxID_slv'],
usecols=['accID','start','end','name','taxID_slv'])

## load NCBI taxIDs w/ species classification
df_ncbimap = pd.read_csv(input_dict['ncbimap'], sep='\t', 
skiprows=1, ## needed for SILVA 138
names=['accID','start','end','path_ncbi','name','taxID_ncbi'],
usecols=['accID','start','end','taxID_ncbi'])

## load NCBI krona database w/ species classification
df_ncbikrona = pd.read_csv(input_dict['ncbikrona'], sep='\t', 
names=['taxID_ncbi','depth_ncbi','targetID_ncbi','rank_ncbi','name'],
usecols=['taxID_ncbi','rank_ncbi'])


## COMPLETE HIGHER TAXA ##

## extactiong additional information
tmp_depth = []
tmp_name = []
tmp_targetID = []
for splt in df_taxlist['pathname_slv'].str.split(';'):
    tmp_depth = tmp_depth + [ len(splt) - 1 ]
    tmp_name = tmp_name + [splt.pop(-2)]
    tmp_targettax = ";".join(splt)
    if tmp_targettax == "":
        tmp_targettax = [1]
    else:
        tmp_targettax = df_taxlist.loc[ df_taxlist['pathname_slv'] == tmp_targettax, 'taxID_slv' ].tolist()
    tmp_targetID = tmp_targetID + tmp_targettax

## assign to df_taxlist
df_taxlist['depth_slv'] = tmp_depth
df_taxlist['name'] = tmp_name
df_taxlist['targetID_slv'] = tmp_targetID

## exclude useless names
df_taxlist.loc[df_taxlist['name'].str.contains("uncultured|unknown|metagenome|unidentified|environmental sample|bacterium enrichment culture|Incertae Sedis|Unknown Family|Unknown Genus"),'rank_slv'] = "no rank"


## COMPLETE SPECIES ##

## merge silva data frames
df_slv = pd.merge(df_slvmap, df_taxlist, how='left', on='taxID_slv')
## rename columns
df_slv.rename(columns = {'pathname_slv':'path_slv', 'name_x':'name'},inplace=True)
## get long accession ID with start and stop
df_slv['accIDstartend'] = df_slv['accID'].map(str) + "." + df_slv['start'].map(str) + "." + df_slv['end'].map(str)

## merge ncbi data frames
df_ncbi = pd.merge(df_ncbimap, df_ncbikrona, how='left', on='taxID_ncbi')
## get long accession ID with start and stop
df_ncbi['accIDstartend'] = df_ncbi['accID'].map(str) + "." + df_ncbi['start'].map(str) + "." + df_ncbi['end'].map(str)
df_ncbi.drop(['accID','start','end'], axis=1, inplace=True)

## merge silva and ncbi data frames
df_slv = pd.merge(df_slv, df_ncbi, how='left', on='accIDstartend')

## new rank for sequences in silva fasta file
df_slv['rank_new'] = "no rank"
## assign rank species to all taxa with genus target
df_slv.loc[ df_slv['rank_slv']=='genus', 'rank_new'] = "species"
## exclude subspecies ranks according to ncbi nomeclature
df_slv.loc[ df_slv['rank_ncbi']=="subspecies", 'rank_new'] = "no rank"
df_slv.loc[ df_slv['rank_ncbi']=="varietas", 'rank_new'] = "no rank"
df_slv.loc[ df_slv['rank_ncbi']=="no rank", 'rank_new'] = "no rank"
## exclude uninformative taxa
df_slv.loc[df_slv['name'].str.contains("uncultured|unknown|metagenome|unidentified|environmental sample|bacterium enrichment culture"),'rank_new'] = "no rank"
df_slv.loc[df_slv['name']=="bacterium",'rank_new'] = "no rank"

## split into used and unused ("no rank") taxa
df_slv_norank = df_slv.loc[ df_slv['rank_new']=="no rank",:]
df_slv = df_slv.loc[ df_slv['rank_new']!="no rank",:]

# ## handle "no rank" targets (higher taxa) for norank species 
# ## create data frame with taxID and rank of target
# df_norank = pd.merge( df_slv_norank.loc[:,['taxID_slv']], df_taxlist, how='left', on='taxID_slv' )
# ## loop until break
# taxID_slv = []
# # targetID_slv = []
# # path_slv = []
# # depth_slv = []
# # name = []
# while True:
#     ## get all taxIDs whoes targetID has "no rank"
#     idx = df_norank['rank_slv']=="no rank"
#     ## if there are non (anymore)
#     if sum(idx) == 0:
#         ## return taxID and targetID
#         taxID_slv = df_norank['taxID_slv']
#         # targetID_slv = df_norank['targetID_slv']
#         # path_slv = df_norank['pathname_slv']
#         # depth_slv = df_norank['depth_slv']
#         # name = []
#         ## break loop
#         break
# 
#     ## for all taxIDs whoes targetID has "no rank", increase the level of the target to the next higher taxon
#     df_norank.loc[idx,'taxID_slv'] = df_norank.loc[idx,'targetID_slv']
#     df_norank = pd.merge( df_norank.loc[:,['taxID_slv']], df_taxlist, how='left', on='taxID_slv' )
# 
# ## assign new values to data frame
# df_slv_norank['taxID_slv'] = taxID_slv
# # df_slv_norank['targetID_slv'] = targetID_slv
# # df_slv_norank['path_slv'] = path_slv
# # df_slv_norank['depth_slv'] = depth_slv
# # df_slv_norank['name'] = name
# del df_norank, taxID_slv #, targetID_slv, path_slv, depth_slv, name

## reconstruct full path with name with tailing ";" at the end
df_slv['pathname_slv'] = df_slv['path_slv'].map(str) + df_slv['name'].map(str) + ";"

## overwrite taxonomy rank for used taxa
df_slv['rank_slv'] = df_slv['rank_new']

## increase depth for used taxa
df_slv['depth_slv'] = df_slv['depth_slv'] + 1

## targetID is taxID for used taxa
df_slv['targetID_slv'] = df_slv['taxID_slv']

## to figure out which accID needs a new taxID, find accIDs with the same depth, rank, target taxID, and name
single_taxIDs = df_slv.loc[:,['depth_slv', 'targetID_slv', 'rank_slv', 'name']].drop_duplicates()
## creata a range of potential taxIDs for the species as long as all higer taxa + new species taxIDs (this should be the maximum number of taxIDs needed, if all would be continuous)
potential_taxIDs = set(range( 10000, len(df_taxlist.index) + len(single_taxIDs.index) + 10000 )) ## +10000 to not use 0 to 9999 as taxID
## exclude all taxIDs used by higer taxa
potential_taxIDs = potential_taxIDs - set(df_taxlist['taxID_slv'])
## reduce to exact number of required taxIDs
potential_taxIDs = list(potential_taxIDs)[0:len(single_taxIDs.index)]
## add new taxIDs to unique species
single_taxIDs['taxID_new'] = potential_taxIDs
## combine with main data frame
df_slv = pd.merge(df_slv, single_taxIDs, how='left', on=['depth_slv', 'targetID_slv', 'rank_slv', 'name'])

## taxID newly created taxID for used taxa
df_slv['taxID_slv'] = df_slv['taxID_new']


## SILVA LIKE ##

## write sequence taxon association
df_tmp = pd.concat([
    ## used taxons with higher taxon ID (was reassignes to targetID here)
    df_slv.loc[:,['accIDstartend','targetID_slv']].rename(columns={"targetID_slv" : "taxID_slv"}),
    ## unused taxons ("no rank") with higher taxon ID (was left as taxID here)
    df_slv_norank.loc[:,['accIDstartend','taxID_slv']]
])
df_tmp.to_csv(
    ## save
    output_dict['krakseq2tax_G'], mode='w', sep='\t', header=False, index=False, quoting=0, float_format="%g"
)
df_tmp.to_csv(
    ## save
    output_dict['kronaseq2tax_G'], mode='w', sep='\t', header=False, index=False, quoting=0, float_format="%g"
)

## write sequence taxon association FOR SPECIES
df_tmp = pd.concat([
    ## used taxons with species taxon ID (was created as taxID here)
    df_slv.loc[:,['accIDstartend','taxID_slv']],
    ## unused taxons ("no rank") with "species" taxon ID (not actually species taxon ID, but left at higer taxID here)
    df_slv_norank.loc[:,['accIDstartend','taxID_slv']]
])
df_tmp.to_csv(
    ## save
    output_dict['krakseq2tax_S'], mode='w', sep='\t', header=False, index=False, quoting=0, float_format="%g"
)
df_tmp.to_csv(
    ## save
    output_dict['kronaseq2tax_S'], mode='w', sep='\t', header=False, index=False, quoting=0, float_format="%g"
)

## remove duplicates: many taxa have more than one accID
df_slv = df_slv.drop_duplicates(subset=['depth_slv', 'targetID_slv', 'rank_slv', 'name'])

## write taxlist-like file (without last columns)
pd.concat([
    ## include root
    pd.DataFrame([[";",1,"root"]],columns=['pathname_slv','taxID_slv','rank_slv']),
    ## include higher ranks
    df_taxlist.loc[:,['pathname_slv','taxID_slv','rank_slv']]
]).to_csv(
    ## save
    output_dict['kronataxlist_G'], mode='w', sep='\t', header=False, index=False, quoting=0, float_format="%g"
)
## write taxlist-like file (without last columns) INCLUDING SPECIES
pd.concat([
    ## include root
    pd.DataFrame([[";",1,"root"]],columns=['pathname_slv','taxID_slv','rank_slv']),
    ## include higher ranks
    df_taxlist.loc[:,['pathname_slv','taxID_slv','rank_slv']],
    ## include species
    df_slv.loc[:,['pathname_slv','taxID_slv','rank_slv']]
]).to_csv(
    ## save
    output_dict['kronataxlist_S'], mode='w', sep='\t', header=False, index=False, quoting=0, float_format="%g"
)


## KRONA ##

## write krona silva reference taxonomy
pd.concat([
    ## include root
    pd.DataFrame([[1,0,1,"no rank","root"]],columns=['taxID_slv','depth_slv','targetID_slv','rank_slv','name']),
    ## include higher ranks
    df_taxlist.loc[:,['taxID_slv','depth_slv','targetID_slv','rank_slv','name']]
]).to_csv(
    ## save
    output_dict['kronataxtab_G'], mode='w', sep='\t', header=False, index=False, quoting=0, float_format="%g"
)
## write krona silva reference taxonomy INCLUDING SPECIES
pd.concat([
    ## include root
    pd.DataFrame([[1,0,1,"no rank","root"]],columns=['taxID_slv','depth_slv','targetID_slv','rank_slv','name']),
    ## include higher ranks
    df_taxlist.loc[:,['taxID_slv','depth_slv','targetID_slv','rank_slv','name']],
    ## include species
    df_slv.loc[:,['taxID_slv','depth_slv','targetID_slv','rank_slv','name']]
]).to_csv(
    ## save
    output_dict['kronataxtab_S'], mode='w', sep='\t', header=False, index=False, quoting=0, float_format="%g"
)


## KRAKEN2 ##

## kraken formats are strange!? Add missing columns
df_taxlist = df_taxlist.assign(
    dummy1 = "|",
    dummy2 = "-",
    dummy3 = "scientific name"
)
df_slv = df_slv.assign(
    dummy1 = "|",
    dummy2 = "-",
    dummy3 = "scientific name"
)

## write kraken2 taxonomy names
pd.concat([
    ## include root
    pd.DataFrame([[1,"|","root","|","-","|","scientific name","|"]],columns=['taxID_slv','dummy1','name','dummy1','dummy2','dummy1','dummy3','dummy1']),
    ## include higher ranks
    df_taxlist.loc[:,['taxID_slv','dummy1','name','dummy1','dummy2','dummy1','dummy3','dummy1']]
]).to_csv(
    ## save
    output_dict['kraknames_G'], mode='w', sep='\t', header=False, index=False, quoting=0, float_format="%g"
)
## write kraken2 taxonomy names INCLUDING SPECIES
pd.concat([
    ## include root
    pd.DataFrame([[1,"|","root","|","-","|","scientific name","|"]],columns=['taxID_slv','dummy1','name','dummy1','dummy2','dummy1','dummy3','dummy1']),
    ## include higher ranks
    df_taxlist.loc[:,['taxID_slv','dummy1','name','dummy1','dummy2','dummy1','dummy3','dummy1']],
    ## include species
    df_slv.loc[:,['taxID_slv','dummy1','name','dummy1','dummy2','dummy1','dummy3','dummy1']]
]).to_csv(
    ## save
    output_dict['kraknames_S'], mode='w', sep='\t', header=False, index=False, quoting=0, float_format="%g"
)

## write kraken2 taxonomy nodes
pd.concat([
    ## include root
    pd.DataFrame([[1,"|",1,"|","no rank","|","-","|"]],columns=['taxID_slv','dummy1','targetID_slv','dummy1','rank_slv','dummy1','dummy2','dummy1']),
    ## include higher ranks
    df_taxlist.loc[:,['taxID_slv','dummy1','targetID_slv','dummy1','rank_slv','dummy1','dummy2','dummy1']]
]).to_csv(
    ## save
    output_dict['kraknodes_G'], mode='w', sep='\t', header=False, index=False, quoting=0, float_format="%g"
)
## write kraken2 taxonomy nodes INCLUDING SPECIES
pd.concat([
    ## include root
    pd.DataFrame([[1,"|",1,"|","no rank","|","-","|"]],columns=['taxID_slv','dummy1','targetID_slv','dummy1','rank_slv','dummy1','dummy2','dummy1']),
    ## include higher ranks
    df_taxlist.loc[:,['taxID_slv','dummy1','targetID_slv','dummy1','rank_slv','dummy1','dummy2','dummy1']],
    ## include species
    df_slv.loc[:,['taxID_slv','dummy1','targetID_slv','dummy1','rank_slv','dummy1','dummy2','dummy1']]
]).to_csv(
    ## save
    output_dict['kraknodes_S'], mode='w', sep='\t', header=False, index=False, quoting=0, float_format="%g"
)
