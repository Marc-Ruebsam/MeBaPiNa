## load library
library(ggplot2)

## load list opf read lengths from fasta index
lenlist <- read.table( snakemake@input[[1]], header = FALSE, sep = "\t",
    col.names = c("names", "length", "V3", "V4", "V5"),
    colClasses = c("NULL","integer","NULL","NULL","NULL"), stringsAsFactors = FALSE )

## calulate occurence
lenfreq <- as.data.frame(table(lenlist))
write.table( lenfreq, file = snakemake@output[["length_stat"]],
    quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE )

## calculate mean, sd, median and mode
lenstats <- c(with( lenlist, c(mean(length), sd(length), median(length))), lenfreq[which.max(lenfreq$Freq),2] )

## plot
lengthplot <- ggplot() + theme_bw() +
    geom_bar( data = lenlist, aes( x = length ) ) +
    geom_vline( aes( xintercept = lenstats[1] ), color = "darkred", size = 0.7 ) +
    geom_polygon( aes( x = c(lenstats[1]-lenstats[2],lenstats[1]-lenstats[2],lenstats[1]+lenstats[2],lenstats[1]+lenstats[2]),
        y = c(-Inf,Inf,Inf,-Inf) ),
        color = "gray40", fill = "gray40", alpha = 0.4 ) +
    geom_label( aes( x = min(lenlist), y = lenstats[4],
        label = paste0("mean=",round(lenstats[1],2),"\nsd=",round(lenstats[2],2),"\nmedian=",round(lenstats[3],2),"\nmode=",round(lenstats[4],2)) ),
        hjust = 0, vjust = 1) +
    xlab("length") + ylab("occurence")

ggsave( snakemake@output[["length_plot"]], lengthplot, units = "cm", height = 12 , width = 21 )
