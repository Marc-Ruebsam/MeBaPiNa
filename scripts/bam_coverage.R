## load library
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(viridisLite))
suppressPackageStartupMessages(library(scales))

## load list of reads per reference position
covlist <- read.table( snakemake@input[[1]], header = FALSE, sep = "\t",
    col.names = c("sequence", "position", "coverage"),
    colClasses = c("NULL","integer","integer"), stringsAsFactors = FALSE )
## assign number to each reference
covlist$seq <- factor(cumsum(covlist$position==1))
## split into list of data frames per reference sequence
covlist <- split(covlist,covlist$seq)
## filter by minumum coverage
covlist <- lapply(covlist,function(covdf){ if( mean(covdf$coverage) >= snakemake@config[["filtering"]][["min_featurereads"]] ){return(covdf)}else{return(NULL)} })
covlist <- covlist[!sapply(covlist,is.null)]

## loop over data frames per reference sequence
n_chunks <- 128
covlist_chunks <- sapply( covlist, function(cov_ref){
  ## calculate mean coverage over all positions
  cov_mean <- mean(cov_ref$coverage)
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

## plot coverage distribution
covdist_plot <- ggplot() + theme_bw() +
  geom_histogram( data = data.frame(coverage=covlist_chunks[[1]]), aes( x = coverage), bins = 100 ) +
  scale_x_log10() +
  xlab("mean reference coverage") + ylab("occurence")
## save
ggsave( snakemake@output[["covdist_plot"]], covdist_plot, units = "cm", height = 12 , width = 21 )

## get median in chunks
cov_mean <- rowMeans(covlist_chunks[[2]])
## melt data frame 
covlist_chunks[[2]]$pos <- seq(0,1,length.out = n_chunks)
covlist_chunks[[2]] <- melt( covlist_chunks[[2]], id.vars = "pos", value.name = "cov" )

covpos_plot <- ggplot() + theme_bw() +
geom_bin2d( data = covlist_chunks[[2]], aes( x = pos, y = cov ), bins = n_chunks ) +
  scale_fill_gradientn( colours = viridis(256, option = "D"), trans = "log", labels = trans_format("identity", function(x) round((x))) ) + 
  geom_line( aes( x = seq(0,1,length.out = n_chunks), y = cov_mean ), color = "red", size = 0.8 ) + 
  xlab("relative position") + ylab("relative coverage")
## save
ggsave( snakemake@output[["covpos_plot"]], covpos_plot, units = "cm", height = 12 , width = 21 )
