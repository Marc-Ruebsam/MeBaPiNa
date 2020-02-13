## SETUP ##

## dependencies
import pandas as pd

## logging
sys.stdout = open(snakemake.log[0], 'w')
sys.stderr = open(snakemake.log[0], 'w')

## input files
input_dict = {
    'taxlist' : snakemake.input['taxlist'],
    'slvmap' : snakemake.input['slvmap'],
    'ncbimap' : snakemake.input['ncbimap'],
    'ncbikrona' : snakemake.input['ncbikrona'] + "/taxonomy.tab"
}

## output files
output_dict = {
    'kraknames_S' : snakemake.output['krak_S'] + "/taxonomy/names.dmp",
    'kraknodes_S' : snakemake.output['krak_S'] + "/taxonomy/nodes.dmp",
    'krakseq2tax_S' : snakemake.output['krak_S'] + "/seqid2taxid.map",
    'kraknames_G' : snakemake.output['krak_G'] + "/taxonomy/names.dmp",
    'kraknodes_G' : snakemake.output['krak_G'] + "/taxonomy/nodes.dmp",
    'krakseq2tax_G' : snakemake.output['krak_G'] + "/seqid2taxid.map",
    'krona_S' : snakemake.output['krona_S'] + "/taxonomy.tab",
    'krona_G' : snakemake.output['krona_G'] + "/taxonomy.tab",
    'taxlist_S' : snakemake.output['taxlist_S'],
    'taxlist_G' : snakemake.output['taxlist_G']
}

# input_dict = {
#     'ncbikrona' : "METADATA/Reference_Sequences/ncbi/krona/taxonomy.tab",
#     'ncbimap' : "METADATA/Reference_Sequences/silva/ncbimap.txt",
#     'slvmap' : "METADATA/Reference_Sequences/silva/slvmap.txt",
#     'taxlist' : "METADATA/Reference_Sequences/silva/taxlist.txt"
# }
# output_dict = {
#     'kraknames_S' : "METADATA/Reference_Sequences/silva/kraken2/species_tmp/taxonomy/names.dmp",
#     'kraknodes_S' : "METADATA/Reference_Sequences/silva/kraken2/species_tmp/taxonomy/nodes.dmp",
#     'krakseq2tax_S' : "METADATA/Reference_Sequences/silva/kraken2/species_tmp/seqid2taxid.map",
#     'kraknames_G' : "METADATA/Reference_Sequences/silva/kraken2/genus_tmp/taxonomy/names.dmp",
#     'kraknodes_G' : "METADATA/Reference_Sequences/silva/kraken2/genus_tmp/taxonomy/nodes.dmp",
#     'krakseq2tax_G' : "METADATA/Reference_Sequences/silva/kraken2/genus_tmp/seqid2taxid.map",
#     'krona_S' : "METADATA/Reference_Sequences/silva/krona/species/taxonomy.tab",
#     'krona_G' : "METADATA/Reference_Sequences/silva/krona/genus/taxonomy.tab",
#     'taxlist_S' : "METADATA/Reference_Sequences/silva/krona/species/taxlist.txt",
#     'taxlist_G' : "METADATA/Reference_Sequences/silva/krona/genus/taxlist.txt"
# }


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
df_taxlist.loc[df_taxlist['name'].isin(["uncultured","Incertae Sedis","Unknown Family"]),'rank_slv'] = "no rank"


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

## reconstruct full path with name with tailing ";" at the end
df_slv['pathname_slv'] = df_slv['path_slv'].map(str) + df_slv['name'].map(str) + ";"

## overwrite taxonomy rank for used taxa
df_slv['rank_slv'] = df_slv['rank_new']

## increase depth for used taxa
df_slv['depth_slv'] = df_slv['depth_slv'] + 1

## targetID is taxID for used taxa
df_slv['targetID_slv'] = df_slv['taxID_slv']

## create custom taxID from silva and ncbi ID for used taxa
df_slv['taxID_new'] = (df_slv['taxID_ncbi'].map(str) + "0" + df_slv['taxID_slv'].map(str)).str.replace('-', '') ## adding additional "0" to wnsure unique ID, also some taxIDs in the ena file are negative
## some taxIDs are used more than once ... -.-
dupli_taxIDs = df_slv.loc[:,['taxID_new', 'depth_slv', 'targetID_slv', 'rank_slv', 'name']].drop_duplicates()['taxID_new'].value_counts()
dupli_taxIDs = dupli_taxIDs[dupli_taxIDs > 1].index.tolist()
## loop over duplicate taxIDs
for dupliID in dupli_taxIDs:
    ## get index of rows matching this id in df_slv
    dupli_rows = df_slv.loc[df_slv['taxID_new']==dupliID,'taxID_new'].index
    ## add a number to the end of the taxID (...0, ...1, ...2)
    unique_taxIDs = [dupliID + str(s) for s in list(range(0,len(dupli_rows)))]
    ## replace duplicate taxID_new with unique taxID_new
    df_slv.loc[dupli_rows,'taxID_new'] = unique_taxIDs

## taxID newly created taxID for used taxa
df_slv['taxID_slv'] = df_slv['taxID_new']


## SILVA LIKE ##

## write sequence taxon association
pd.concat([
    ## used taxons with higher taxon ID (was reassignes to targetID here)
    df_slv.loc[:,['accIDstartend','targetID_slv']],
    ## unused taxons ("no rank") with higher taxon ID (was left as taxID here)
    df_slv_norank.loc[:,['accIDstartend','taxID_slv']],
]).to_csv(
    ## save
    output_dict['krakseq2tax_G'], mode='w', sep='\t', header=False, index=False, quoting=0, float_format="%g"
)
## write sequence taxon association FOR SPECIES
pd.concat([
    ## used taxons with species taxon ID (was created as taxID here)
    df_slv.loc[:,['accIDstartend','taxID_slv']],
    ## unused taxons ("no rank") with "species" taxon ID (not actually species taxon ID, but left at higer taxID here)
    df_slv_norank.loc[:,['accIDstartend','taxID_slv']],
]).to_csv(
    ## save
    output_dict['krakseq2tax_S'], mode='w', sep='\t', header=False, index=False, quoting=0, float_format="%g"
)

## remove duplicates: many taxa have more than one accID
df_slv = df_slv.drop_duplicates(subset=['taxID_slv','depth_slv','targetID_slv','rank_slv','name'])

## write taxlist-like file (without last columns)
pd.concat([
    ## include root
    pd.DataFrame([[1,";","root"]],columns=['taxID_slv','pathname_slv','rank_slv']),
    ## include higher ranks
    df_taxlist.loc[:,['taxID_slv','pathname_slv','rank_slv']]
]).to_csv(
    ## save
    output_dict['taxlist_G'], mode='w', sep='\t', header=False, index=False, quoting=0, float_format="%g"
)
## write taxlist-like file (without last columns) INCLUDING SPECIES
pd.concat([
    ## include root
    pd.DataFrame([[1,";","root"]],columns=['taxID_slv','pathname_slv','rank_slv']),
    ## include higher ranks
    df_taxlist.loc[:,['taxID_slv','pathname_slv','rank_slv']],
    ## include species
    df_slv.loc[:,['taxID_slv','pathname_slv','rank_slv']]
]).to_csv(
    ## save
    output_dict['taxlist_S'], mode='w', sep='\t', header=False, index=False, quoting=0, float_format="%g"
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
    output_dict['krona_G'], mode='w', sep='\t', header=False, index=False, quoting=0, float_format="%g"
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
    output_dict['krona_S'], mode='w', sep='\t', header=False, index=False, quoting=0, float_format="%g"
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
