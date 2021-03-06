## PACKAGES ##

## load packages
suppressPackageStartupMessages(library("phyloseq"))
suppressPackageStartupMessages(library("vegan"))
suppressPackageStartupMessages(library("ggplot2"))

## set variables
below_flt <- snakemake@config[["filtering"]][["min_featurereads"]]
above_flt <- 10000 ## group references above or equal to this mean (rounded down) coverage

## LOAD INPUT ##

# ## load meta data
# meta_df <- read.csv( snakemake@input[["taxlist"]], sep="\t", header = FALSE, fill = TRUE, stringsAsFactors = FALSE,
#     ## we know the types of content to expect
#     colClasses = c( "character", "integer", "factor" ),
#     ## create column names
#     col.names = c( "path", "taxid", "rank" ) )

## get number of rows and columns in sample file
cols <- as.numeric(system(paste0( "awk -F '\\t' 'BEGIN{numb=0; OFS = \"\\n\"}; {if (NF>numb){numb=NF}}; END{print numb,NR}' ", snakemake@input[["sample_file"]] ), intern = TRUE))
## table might be empty
if( cols[2] == 0 ){
    div_df <- matrix(c("NA","richness_S",
        "NA","shannon_H_base_e",
        "NA","shannon_H_base_2",
        "NA","simpson_lambda",
        "NA","simpson_lambda_inverse",
        "NA","simpson_lambda_rarefy",
        "NA","evenness_J",
        "NA","fishers_alpha"),ncol=2,byrow=TRUE)
    ## return
    write.table(div_df, file = snakemake@output[["report"]], quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

    ## plot
    covdist_plot <- ggplot() + theme_bw() +
      geom_label( aes( x = above_flt, y = 1, label = "no data passed thresholds" ), hjust = 1, vjust = 1) +
      xlab("taxa abundance") + ylab("occurence") +
      scale_x_log10()
    ## save
    ggsave( snakemake@output[["covdist_plot"]], covdist_plot, units = "cm", height = 12 , width = 21 )

    ## exit
    quit(save = "no", status = 0)
}
## import table
input_df <- read.csv(snakemake@input[["sample_file"]], sep="\t", header = FALSE, fill = TRUE, stringsAsFactors = FALSE,
    ## we know the types of content to expect
    colClasses = c( "integer", rep("character", (cols[1]-1)) ),
    ## create column names
    col.names = c( "counts", paste0("depth", c(1:(cols[1]-1))) ),
    ## we know the number of rows to expect
    nrows = cols[2] )
## curration
input_df[input_df == ""] <- NA
input_df <- input_df[input_df$counts != 0,]

## name of last rank should all be species
specs_name <- apply( input_df, 1, function(row){ row[max(which(!is.na(row)))] } )
attr(specs_name,"names") <- NULL
# ## split counts and taxonomic paths (path list)
# counts <- input_df[,"counts"]
# specs_path_lst <- input_df[,which(!colnames(input_df) %in% "counts")]
# ## collapse path to single string
# specs_path <- paste0(apply( specs_path_lst, 1, function(specs){
#     specs <- specs[!is.na(specs)]
#     paste(specs,  sep = "", collapse = ";")
# } ),";")
# ## reverse order of ranks in path list
# specs_path_lst <- specs_path_lst[,rev(colnames(specs_path_lst))]

## count table
counts_df <- data.frame(t(input_df$counts))
colnames(counts_df) <- specs_name


## PLOT ##

## get count values
plot_df <- input_df$counts
## get numbers max coverage number
plot_flt <- c( max = max(plot_df) )
## group values above threshold
idx <- plot_df >= above_flt
plot_flt["above"] <- sum(idx)
plot_df[ idx ] <- above_flt

## plot coverage distribution
if( length(plot_df) != 0 ){
  covdist_plot <- ggplot() + theme_bw() +
    geom_bar( data = data.frame(abundance=plot_df), aes( x = abundance), width = 0.02 ) +
    geom_label( aes( x = max(plot_df), y = max(table(plot_df)),
      label = paste0("exclude: taxa < ",below_flt,"\ngroup: ",plot_flt["above"]," taxa >= ",above_flt,"\nwith max cov of ",plot_flt["max"]) ),
      hjust = 1, vjust = 1) +
    xlab("taxa abundance") + ylab("occurence") +
    scale_x_log10()
}else{
  covdist_plot <- ggplot() + theme_bw() +
    geom_label( aes( x = above_flt, y = 1, label = "no data passed thresholds" ), hjust = 1, vjust = 1) +
    xlab("taxa abundance") + ylab("occurence") +
    scale_x_log10()
}
## save
ggsave( snakemake@output[["covdist_plot"]], covdist_plot, units = "cm", height = 12 , width = 21 )


## DIVERSITY INDICES ##

## get common diversity indices
div_df <- rbind(
    ## pecies richness (S)
    c(as.character(
        specnumber(counts_df) ),
        "richness_S"),
    ## Shannon–Weaver
    ## The Shannon entropy quantifies the uncertainty in predicting the species identity of an individual that is taken at random from the dataset.
    c(as.character(round(
        diversity(counts_df, index = "shannon", base = exp(1)) ,2)),
        "shannon_H_base_e"),
    ## shannon base 2
    c(as.character(round(
        diversity(counts_df, index = "shannon", base = 2) ,2)),
        "shannon_H_base_2"),
    ## Simpson
    ## The measure equals the probability that two entities taken at random from the dataset of interest represent the same type
    c(as.character(round(
        diversity(counts_df, index = "simpson") ,2)),
        "simpson_lambda"),
    ## inverse Simpson
    c(as.character(round(
        diversity(counts_df, index = "invsimpson") ,2)),
        "simpson_lambda_inverse"),
    ## Unbiased Simpson (Hurlbert 1971, eq. 5) with rarefy
    #!# rarefying is downsampling. should not be used when proper normalization is applied
    c(as.character(round(
        (rarefy(counts_df, 2) - 1) ,2)),
        "simpson_lambda_rarefy"),
    ## Pielou's evenness (J):
    c(as.character(round(
        diversity(counts_df, index = "shannon", base = exp(1))/log(specnumber(counts_df)) ,2)),
        "evenness_J"),
    ## alpha parameter of Fisher’s log-series
    c(as.character(round(
        fisher.alpha(counts_df) ,2)),
        "fishers_alpha")
)
## save
write.table(div_df, file = snakemake@output[["report"]], quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


# ## TAXONOMIC DIVERSITY ##
#
# null <- sapply( dimnames(sample_lst)[[2]], function(sample_name){
#     name_path_lst <- sample_lst[,sample_name][[3]]
#     name_path_lst[,2:ncol(name_path_lst)]
#     name_taxdist <- taxa2dist( name_path_lst, varstep = TRUE )
#     taxondive(dune, taxdis)
#     str(name_taxdist)
# })
# taxa2dist(sample_lst[3,][[6]][,2:], varstep=TRUE)
#
# ## taxonomic diversity, taxonomic distinctness, ...
# taxdis <- taxa2dist(dune.taxon, varstep=TRUE)
# tr <- hclust(taxdis, "aver")
# ## Taxonomic diversity: average distance of trai
# taxondive(dune, taxdis)
# ## Functional diversity: the height of trait tree
# treedive(dune, tr)
#
# ## SPECIES ABUNDANCE MODELS ##
#
# fisherfit(BCI[1,])
# prestondistr(BCI[1,])
# radfit(BCI[1,])
#
# ## SPECIES ACCUMULATION AND BETA DIVERSITY ##
#
# ## beta diversity defined as gamma/alpha - 1:
# #!# between all samples
# alpha <- with(dune.env, tapply(specnumber(dune), Management, mean))
# gamma <- with(dune.env, specnumber(dune, Management))
# gamma/alpha - 1
# ncol(BCI)/mean(specnumber(BCI)) - 1
#
# ## Sorensen index of dissimilarity (pairwise)
# beta <- vegdist(BCI, binary=TRUE)
# mean(beta)
#
# ## all sorts of other methods
# betadiver(BCI)
