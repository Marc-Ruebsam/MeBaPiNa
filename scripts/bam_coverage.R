## load library
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(viridisLite))
suppressPackageStartupMessages(library(scales))

## set variables
n_chunks <- 128 ## number of chunks to cut reference sequences into
below_flt <- snakemake@config[["filtering"]][["min_featurereads"]] ## exclude reference below this mean (rounded down) coverage
above_flt <- 1000 ## group references above or equal to this mean (rounded down) coverage

## load list of reads per reference position
covlist <- read.table( snakemake@input[[1]], header = FALSE, sep = "\t",
    col.names = c("sequence", "position", "coverage"),
    colClasses = c("NULL","integer","integer"), stringsAsFactors = FALSE )
## assign number to each reference
covlist$seq <- factor(cumsum(covlist$position==1))
## split into list of data frames per reference sequence
covlist <- split(covlist,covlist$seq)

## loop over data frames per reference sequence
covlist_chunks <- sapply( covlist, function(cov_ref){
  ## calculate mean coverage over all positions (rounded down)
  cov_mean <- floor(mean(cov_ref$coverage))
  ## mean coverage in chunks
  cov_mean_chunk  <- sapply( split( cov_ref$coverage, cut(cov_ref$position,n_chunks) ), mean )
  ## normalize to max coverage
  cov_mean_chunk <- cov_mean_chunk/max(cov_mean_chunk)
  ## return
  return(list(cov_mean,cov_mean_chunk))
} )
## restructure
covlist_chunks <- list( unlist(covlist_chunks[1,]),
  data.frame(covlist_chunks[2,]) )

## get numbers max coverage number
cov_flt <- c( max = max(covlist_chunks[[1]]) )
## filter references below threshold
idx <- covlist_chunks[[1]] < below_flt
cov_flt["below"] <- sum(idx)
covlist_chunks[[1]] <- covlist_chunks[[1]][!idx]
covlist_chunks[[2]] <- covlist_chunks[[2]][!idx] ## also for chunks
## group values above threshold
idx <- covlist_chunks[[1]] >= above_flt
cov_flt["above"] <- sum(idx)
covlist_chunks[[1]][ idx ] <- 1000

## get median in chunks
cov_mean <- rowMeans(covlist_chunks[[2]])
## melt data frame
covlist_chunks[[2]]$pos <- seq(0,1,length.out = n_chunks)
covlist_chunks[[2]] <- melt( covlist_chunks[[2]], id.vars = "pos", value.name = "cov" )

## plot coverage distribution
if( length(covlist_chunks[[1]]) != 0 ){
  covdist_plot <- ggplot() + theme_bw() +
    geom_bar( data = data.frame(coverage=covlist_chunks[[1]]), aes( x = coverage), width = 0.02 ) +
    geom_label( aes( x = max(covlist_chunks[[1]]), y = max(table(covlist_chunks[[1]])),
      label = paste0("exclude: ",cov_flt["below"]," refs < ",below_flt,"\ngroup: ",cov_flt["above"]," refs >= ",above_flt,"\nwith max cov of ",cov_flt["max"]) ),
      hjust = 1, vjust = 1) +
    xlab("mean reference coverage") + ylab("occurence") +
    scale_x_log10()
}else{
  covdist_plot <- ggplot() + theme_bw() +
    geom_label( aes( x = above_flt, y = 1, label = "no data passed thresholds" ), hjust = 1, vjust = 1) +
    xlab("mean reference coverage") + ylab("occurence") +
    scale_x_log10()
}
## save
Sys.sleep(3) #!# attempt to correctly save plots
ggsave( snakemake@output[["covdist_plot"]], covdist_plot, units = "cm", height = 12 , width = 21 )
Sys.sleep(3) #!# attempt to correctly save plots

## plot positional coverage
if( length(covlist_chunks[[1]]) != 0 ){
  covpos_plot <- ggplot() + theme_bw() +
    geom_bin2d( data = covlist_chunks[[2]], aes( x = pos, y = cov ), bins = n_chunks ) +
    scale_fill_gradientn( colours = viridis(256, option = "D"), trans = "log", labels = trans_format("identity", function(x) round((x))) ) +
    geom_line( aes( x = seq(0,1,length.out = n_chunks), y = cov_mean ), color = "darkred", size = 0.7 ) +
    xlab("relative position") + ylab("relative coverage")
}else{
  covpos_plot <- ggplot() + theme_bw() +
    geom_label( aes( x = 0.5, y = 1, label = "no data passed thresholds" ), hjust = 1, vjust = 1) +
    xlab("relative position") + ylab("relative coverage")
}
## save
ggsave( snakemake@output[["covpos_plot"]], covpos_plot, units = "cm", height = 12 , width = 21 )
