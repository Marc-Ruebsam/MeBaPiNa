## INPUT ##

library("phyloseq")
library("DESeq2")
library("vegan")


## LOAD INPUT ##

## meta data
meta_files <- c(
    taxlist = paste0( "METADATA/Reference_Sequences/silva/krona/species/", "taxlist.txt" ),
    metadata = "METADATA/EXPERIMENT_SEQUENCING.xlsx"
)

## sample data
sample_files <- c(
    "01_processed_data/03_kmer_mapping/20191007_1559_MN31344_FAK76605_2bf006ff/barcode01/silva_species/filtered.ktaxlist",
    "01_processed_data/03_kmer_mapping/20191007_1559_MN31344_FAK76605_2bf006ff/barcode02/silva_species/filtered.ktaxlist",
    "01_processed_data/03_kmer_mapping/20191007_1559_MN31344_FAK76605_2bf006ff/barcode03/silva_species/filtered.ktaxlist",
    "01_processed_data/03_kmer_mapping/20191007_1559_MN31344_FAK76605_2bf006ff/barcode04/silva_species/filtered.ktaxlist",
    "01_processed_data/03_kmer_mapping/20191007_1559_MN31344_FAK76605_2bf006ff/barcode05/silva_species/filtered.ktaxlist",
    "01_processed_data/03_kmer_mapping/20191007_1559_MN31344_FAK76605_2bf006ff/barcode06/silva_species/filtered.ktaxlist",
    "01_processed_data/03_kmer_mapping/20191007_1559_MN31344_FAK76605_2bf006ff/barcode07/silva_species/filtered.ktaxlist",
    "01_processed_data/03_kmer_mapping/20191007_1559_MN31344_FAK76605_2bf006ff/barcode08/silva_species/filtered.ktaxlist",
    "01_processed_data/03_kmer_mapping/20191007_1559_MN31344_FAK76605_2bf006ff/barcode09/silva_species/filtered.ktaxlist",
    "01_processed_data/03_kmer_mapping/20191007_1559_MN31344_FAK76605_2bf006ff/barcode10/silva_species/filtered.ktaxlist",
    "01_processed_data/03_kmer_mapping/20191007_1559_MN31344_FAK76605_2bf006ff/barcode11/silva_species/filtered.ktaxlist",
    "01_processed_data/03_kmer_mapping/20191007_1559_MN31344_FAK76605_2bf006ff/barcode12/silva_species/filtered.ktaxlist"
)

## load meta data
meta_df <- read.csv( meta_files["taxlist"], sep="\t", header = FALSE, fill = TRUE, stringsAsFactors = FALSE,
    ## we know the types of content to expect
    colClasses = c( "character", "integer", "factor" ),
    ## create column names
    col.names = c( "path", "taxid", "rank" ) )

## load all sample data
sample_lst <- sapply(sample_files, function(input){
    ## get number of rows and columns in file
    cols <- as.numeric(system(paste0( "awk -F '\\t' 'BEGIN{numb=0; OFS = \"\\n\"}; {if (NF>numb){numb=NF}}; END{print numb,NR}' ", input ), intern = TRUE))
    ## table might be empty
    if( cols[2] == 0 ){ 
        return( list( character(0), numeric(0), data.frame(
        depth1=character(0),depth2=character(0),depth3=character(0),depth4=character(0),depth5=character(0),depth6=character(0),depth7=character(0),
        stringsAsFactors=FALSE ), character(0) ) )
    }
    ## import table
    input_df <- read.csv(input, sep="\t", header = FALSE, fill = TRUE, stringsAsFactors = FALSE,
        ## we know the types of content to expect
        colClasses = c( "integer", rep("character", (cols[1]-1)) ),
        ## create column names
        col.names = c( "counts", paste0("depth", c(1:(cols[1]-1))) ),
        ## we know the number of rows to expect
        nrows = cols[2] )
    input_df[input_df == ""] <- NA
    input_df <- input_df[input_df$counts != 0,]

    specs_name <- apply( input_df, 1, function(row){ row[max(which(!is.na(row)))] } )
    attr(specs_name,"names") <- NULL
    counts <- input_df[,"counts"]
    specs_path_lst <- input_df[,which(!colnames(input_df) %in% "counts")]
    specs_path <- paste0(apply( specs_path_lst, 1, function(specs){
        specs <- specs[!is.na(specs)]
        paste(specs,  sep = "", collapse = ";")
    } ),";")
    specs_path_lst <- specs_path_lst[,rev(colnames(specs_path_lst))]

    list( specs_name, counts, specs_path_lst, specs_path )
})
dimnames(sample_lst) <- list(NULL,sapply(strsplit(dimnames(sample_lst)[[2]],"/"),"[",4)) #!#

## initialize count table for all species and all samples
counts_df <- as.data.frame(matrix( data = 0, nrow = ncol(sample_lst), ncol = nrow(meta_df), dimnames = list( dimnames(sample_lst)[[2]], meta_df[,"path"] ) ),stringsAsFactors=FALSE)
## fill count table with read counts
null <- sapply( dimnames(sample_lst)[[2]], function(sample_name){
    name_lst <- sample_lst[,sample_name]
    if( length(name_lst[[1]]) == 0 ){ return(sample_name) }
    counts_df[sample_name, match( name_lst[[4]], colnames(counts_df) )] <<- name_lst[[2]]
    sample_name
} )
## remove species not present in the samples
idx <- colSums(counts_df) != 0
counts_df <- counts_df[,idx]
meta_df <- meta_df[idx,]

## DIVERSITY INDICES ##

## get common diversity indices
## pecies richness (S)
sample_div <- data.frame(richness_S = specnumber(counts_df))

## Shannon–Weaver
sample_div[,"shannon_H"] <- diversity(counts_df, index = "shannon", base = exp(1))
sample_div[,"shannon_base2"] <- diversity(counts_df, index = "shannon", base = 2)
## Simpson
sample_div[,"simp"] <- diversity(counts_df, index = "simpson")
## inverse Simpson
sample_div[,"simp_iverse"] <- diversity(counts_df, index = "invsimpson")
## Unbiased Simpson (Hurlbert 1971, eq. 5) with rarefy
#!# rarefying is downsampling. should not be used when proper normalization is applied
sample_div[,"simp_unbias"] <- rarefy(counts_df, 2) - 1

## Pielou's evenness (J):
sample_div[,"evenness_J"] <- sample_div[,"shannon_H"]/log(sample_div[,"richness_S"])

## alpha parameter of Fisher’s log-series
sample_div[,"alpha"] <- fisher.alpha(counts_df)


## TAXONOMIC DIVERSITY ##

















NEEDS ALL TAXA


null <- sapply( dimnames(sample_lst)[[2]], function(sample_name){
    name_path_lst <- sample_lst[,sample_name][[3]]
    name_path_lst[,2:ncol(name_path_lst)]
    name_taxdist <- taxa2dist( name_path_lst, varstep = TRUE )
    taxondive(dune, taxdis)
    str(name_taxdist)
})
taxa2dist(sample_lst[3,][[6]][,2:], varstep=TRUE)




## taxonomic diversity, taxonomic distinctness, ...
taxdis <- taxa2dist(dune.taxon, varstep=TRUE)
tr <- hclust(taxdis, "aver")
## Taxonomic diversity: average distance of trai
taxondive(dune, taxdis)
## Functional diversity: the height of trait tree
treedive(dune, tr)


## SPECIES ABUNDANCE MODELS ##

fisherfit(BCI[1,])
prestondistr(BCI[1,])
radfit(BCI[1,])


## SPECIES ACCUMULATION AND BETA DIVERSITY ##

## beta diversity defined as gamma/alpha - 1:
#!# between all samples
alpha <- with(dune.env, tapply(specnumber(dune), Management, mean))
gamma <- with(dune.env, specnumber(dune, Management))
gamma/alpha - 1
ncol(BCI)/mean(specnumber(BCI)) - 1

## Sorensen index of dissimilarity (pairwise)
beta <- vegdist(BCI, binary=TRUE)
mean(beta)

## all sorts of other methods
betadiver(BCI)
