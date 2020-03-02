## SETUP ##

## dependencies
import pandas as pd

## logging
sys.stdout = open(snakemake.log[0], 'w')
sys.stderr = open(snakemake.log[0], 'w')

## input files
input_dict = {
    'taxlist' : snakemake.input['taxlist'],
    'slvmap' : snakemake.input['slvmap']
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
    'kronaseq2tax_G' : snakemake.output['kronaseq2tax_G'],
    'qiimetax_S' : snakemake.output['qiimetax_S'],
    'qiimetax_G' : snakemake.output['qiimetax_G']
}

# input_dict = {
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
#     'qiimetax_G' : "METADATA/Reference_Sequences/silva/qiime/genus/taxonomy.tsv",
#     'qiimetax_S' : "METADATA/Reference_Sequences/silva/qiime/species/taxonomy.tsv"
# }
# import os
# krak_dir_S = "/".join(output_dict['krakseq2tax_S'].split("/")[:-1])
# krak_dir_G = "/".join(output_dict['krakseq2tax_G'].split("/")[:-1])
# krona_dir_S = "/".join(output_dict['kronataxtab_S'].split("/")[:-1])
# krona_dir_G = "/".join(output_dict['kronataxtab_G'].split("/")[:-1])
# qiime_dir_S = "/".join(output_dict['qiimetax_S'].split("/")[:-1])
# qiime_dir_G = "/".join(output_dict['qiimetax_G'].split("/")[:-1])
# os.makedirs( krak_dir_S + "/library" )
# os.makedirs( krak_dir_S + "/taxonomy" )
# os.makedirs( krak_dir_G + "/library" )
# os.makedirs( krak_dir_G + "/taxonomy" )
# os.makedirs( krona_dir_S )
# os.makedirs( krona_dir_G )
# os.makedirs( qiime_dir_S )
# os.makedirs( qiime_dir_G )


## LOAD DATA ##

## load taxID list
df_taxlist = pd.read_csv(input_dict['taxlist'], sep='\t', 
names=['pathname','taxID','rank','remark','release'], 
usecols=['pathname','taxID','rank'])

## load SILVA taxIDs w/o species classification
df_accmap = pd.read_csv(input_dict['slvmap'], sep='\t', 
skiprows=1,
names=['accID','start','end','path','name','taxID'],
usecols=['accID','start','end','name','taxID'])


## PER ACCID ##

## add path to accmap
df_accmap = pd.merge(df_accmap, df_taxlist.loc[:,['pathname','taxID']], how='left', on='taxID').rename(columns={'pathname' : 'path', 'taxID' : 'taxID_old'})
## get long accession ID with start and stop
df_accmap['accIDstartend'] = df_accmap['accID'].map(str) + "." + df_accmap['start'].map(str) + "." + df_accmap['end'].map(str)
df_accmap.drop(['accID','start','end'], axis=1, inplace=True)
## set accIDstartend as index
df_accmap.set_index('accIDstartend',inplace=True)
## replace any unwated characters in species name and split the name at widespace
name_fix = df_accmap['name'].str.replace('[^\w\-\[\]\(\)\\/. ]', '').str.split()
#!# concatenate the first two words (or less) of species name back together
df_accmap['name'] = name_fix.str[0].str.cat( name_fix.str[1], sep=" ", na_rep='' ).str.strip()

## get data frame with names from path
df_pathname = df_accmap['path'].str.split(';', expand=True)
#!# remove any unwanted columns
df_pathname.drop(columns=range(7,len(df_pathname.columns)),inplace=True)
## add species name to data frame of names
df_pathname[6] = df_accmap['name']
## replace None with empty string
df_pathname.fillna("", inplace=True)
## add dummy string to other missing taxa
df_pathname.loc[df_pathname[1]=="",1] = "unknown_phylum"
df_pathname.loc[df_pathname[2]=="",2] = "unknown_class"
df_pathname.loc[df_pathname[3]=="",3] = "unknown_order"
df_pathname.loc[df_pathname[4]=="",4] = "unknown_family"
df_pathname.loc[df_pathname[5]=="",5] = "unknown_genus"
## add semicolon to names for path formation
df_pathname = df_pathname + ";"
## combine names to paths
for n in range(1,len(df_pathname.columns)):
    df_pathname[n] = df_pathname[n-1] + df_pathname[n]

## set accIDstartend back to column of accmap
df_accmap.reset_index(inplace=True)


## PER TAXID ##

## set pathname to index of taxlist
df_taxlist.set_index('pathname',inplace=True)

## loop over ranks to:
##     create new taxIDs for new taxa
##     add new taxa to taxlist (including species)
for n in range(0,len(df_pathname.columns)):
    ## get unknown taxa from path (should only be the "unknown_" ones  we just created and all species taxa)
    missing_pathname = set(df_pathname[n]) - set(df_taxlist.index)
    ## skip this rank if no missing taxa are detected
    if missing_pathname == set():
        continue
    ## create new taxIDs above the current max taxID
    new_taxID = list(range( df_taxlist['taxID'].max() + 1, df_taxlist['taxID'].max() + len(missing_pathname) + 1 ))
    ## rank will be set to no rank, except for species
    current_rank = "species" if (n == 6) else "no rank"
    ## add new taxa to taxlist
    df_taxlist = df_taxlist.append( pd.DataFrame({ 'taxID' : new_taxID, 'rank' : current_rank }, index = pd.Index( missing_pathname, name = 'pathname' )) )

## add target taxIDs to taxlist (temporarily sets targetID of domains to their taxID, see next step)
df_taxlist['targetID'] = df_taxlist.loc[df_taxlist.index.str.rsplit(";",n=2).str[0] + ";",'taxID'].to_list()
## set target of domains to root (taxID 1)
df_taxlist.loc[ df_taxlist['taxID'] == df_taxlist['targetID'], 'targetID' ] = 1

## add new taxIDs to accmap
df_accmap['taxID'] = df_taxlist.loc[df_pathname[6],'taxID'].to_list()

## set pathname back to column of taxlist
df_taxlist.reset_index(inplace=True)

## add depth to taxlist
df_taxlist['depth'] = df_taxlist['pathname'].str.split(";").str.len()-1
## add name to taxlist
df_taxlist['name'] = df_taxlist['pathname'].str.split(";").str[-2]


## EXPORT ACCID-TAXID ASSOCIATION ##

## write sequence taxon association FOR HIGHER TAXA
df_tmp = df_accmap.loc[:,['accIDstartend','taxID_old']] ## use original association to higher taxa
df_tmp.to_csv(
    ## save
    output_dict['krakseq2tax_G'], mode='w', sep='\t', header=False, index=False, quoting=0, float_format="%g"
)
df_tmp.to_csv(
    ## save
    output_dict['kronaseq2tax_G'], mode='w', sep='\t', header=False, index=False, quoting=0, float_format="%g"
)
del df_tmp

## write sequence taxon association INCLUDING SPECIES
df_tmp = df_accmap.loc[:,['accIDstartend','taxID']] ## use new taxIDs specific for species
df_tmp.to_csv(
    ## save
    output_dict['krakseq2tax_S'], mode='w', sep='\t', header=False, index=False, quoting=0, float_format="%g"
)
df_tmp.to_csv(
    ## save
    output_dict['kronaseq2tax_S'], mode='w', sep='\t', header=False, index=False, quoting=0, float_format="%g"
)
del df_tmp


## EXPORT TAXLIST ##

## write taxlist-like file (without last columns) FOR HIGHER TAXA
pd.concat([
    ## include root
    pd.DataFrame([[";",1,"root"]],columns=['pathname','taxID','rank']),
    ## include higher ranks
    df_taxlist.loc[ df_taxlist['rank']!="species", ['pathname','taxID','rank'] ] ## exclude species
]).to_csv(
    ## save
    output_dict['kronataxlist_G'], mode='w', sep='\t', header=False, index=False, quoting=0, float_format="%g"
)
## write taxlist-like file (without last columns) INCLUDING SPECIES
pd.concat([
    ## include root
    pd.DataFrame([[";",1,"root"]],columns=['pathname','taxID','rank']),
    ## include higher ranks and species
    df_taxlist.loc[:,['pathname','taxID','rank']]
]).to_csv(
    ## save
    output_dict['kronataxlist_S'], mode='w', sep='\t', header=False, index=False, quoting=0, float_format="%g"
)


## EXPORT KRONA ##

## write krona silva reference taxonomy FOR HIGHER TAXA
pd.concat([
    ## include root
    pd.DataFrame([[1,0,1,"no rank","root"]],columns=['taxID','depth','targetID','rank','name']),
    ## include higher ranks
    df_taxlist.loc[ df_taxlist['rank']!="species", ['taxID','depth','targetID','rank','name'] ]
]).to_csv(
    ## save
    output_dict['kronataxtab_G'], mode='w', sep='\t', header=False, index=False, quoting=0, float_format="%g"
)
## write krona silva reference taxonomy INCLUDING SPECIES
pd.concat([
    ## include root
    pd.DataFrame([[1,0,1,"no rank","root"]],columns=['taxID','depth','targetID','rank','name']),
    ## include higher ranks and species
    df_taxlist.loc[:,['taxID','depth','targetID','rank','name']]
]).to_csv(
    ## save
    output_dict['kronataxtab_S'], mode='w', sep='\t', header=False, index=False, quoting=0, float_format="%g"
)


## EXPORT KRAKEN2 ##

## kraken formats are strange!? Add missing columns
df_taxlist = df_taxlist.assign(
    dummy1 = "|",
    dummy2 = "-",
    dummy3 = "scientific name"
)

## write kraken2 taxonomy names FOR HIGHER TAXA
pd.concat([
    ## include root
    pd.DataFrame([[1,"|","root","|","-","|","scientific name","|"]],columns=['taxID','dummy1','name','dummy1','dummy2','dummy1','dummy3','dummy1']),
    ## include higher ranks
    df_taxlist.loc[ df_taxlist['rank']!="species", ['taxID','dummy1','name','dummy1','dummy2','dummy1','dummy3','dummy1'] ]
]).to_csv(
    ## save
    output_dict['kraknames_G'], mode='w', sep='\t', header=False, index=False, quoting=0, float_format="%g"
)
## write kraken2 taxonomy names INCLUDING SPECIES
pd.concat([
    ## include root
    pd.DataFrame([[1,"|","root","|","-","|","scientific name","|"]],columns=['taxID','dummy1','name','dummy1','dummy2','dummy1','dummy3','dummy1']),
    ## include higher ranks and species
    df_taxlist.loc[:,['taxID','dummy1','name','dummy1','dummy2','dummy1','dummy3','dummy1']]
]).to_csv(
    ## save
    output_dict['kraknames_S'], mode='w', sep='\t', header=False, index=False, quoting=0, float_format="%g"
)

## write kraken2 taxonomy nodes FOR HIGHER TAXA
pd.concat([
    ## include root
    pd.DataFrame([[1,"|",1,"|","no rank","|","-","|"]],columns=['taxID','dummy1','targetID','dummy1','rank','dummy1','dummy2','dummy1']),
    ## include higher ranks
    df_taxlist.loc[ df_taxlist['rank']!="species", ['taxID','dummy1','targetID','dummy1','rank','dummy1','dummy2','dummy1'] ]
]).to_csv(
    ## save
    output_dict['kraknodes_G'], mode='w', sep='\t', header=False, index=False, quoting=0, float_format="%g"
)
## write kraken2 taxonomy nodes INCLUDING SPECIES
pd.concat([
    ## include root
    pd.DataFrame([[1,"|",1,"|","no rank","|","-","|"]],columns=['taxID','dummy1','targetID','dummy1','rank','dummy1','dummy2','dummy1']),
    ## include higher ranks and species
    df_taxlist.loc[:,['taxID','dummy1','targetID','dummy1','rank','dummy1','dummy2','dummy1']]
]).to_csv(
    ## save
    output_dict['kraknodes_S'], mode='w', sep='\t', header=False, index=False, quoting=0, float_format="%g"
)


## EXPORT QIIME ##

## write qiime accID to path FOR HIGHER TAXA
df_tmp = df_pathname[6].str.replace("\s","_").str.split(";") 
df_save = (
    "d__" + df_tmp.str[0] + 
    "; p__" + df_tmp.str[1] + 
    "; c__" + df_tmp.str[2] + 
    "; o__" + df_tmp.str[3] + 
    "; f__" + df_tmp.str[4] + 
    "; g__" + df_tmp.str[5]
)
df_save.reset_index().to_csv(
    ## save
    output_dict['qiimetax_G'], mode='w', sep='\t', header=False, index=False, quoting=0, float_format="%g"
)
## write qiime accID to path INCLUDING SPECIES
(
    df_save + 
    "; s__" + df_tmp.str[6]
).reset_index().to_csv(
    ## save
    output_dict['qiimetax_S'], mode='w', sep='\t', header=False, index=False, quoting=0, float_format="%g"
)
