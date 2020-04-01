library(ggplot2)
library(reshape2)
library(RColorBrewer)

ref_demux_filt_counts <- c(44, 30, 673100, 1487303, 718362, 857069, 678881, 530264, 1341020, 988081, 724721, 625320, 66260, 283913, 723653, 191520, 192578, 157178, 138630, 61976, 78003, 158098, 230634, 932741)

# ## k-mer
# I_03 <- list(
#  total = 672290,
#  species = c(
#   301,
#   22360,
#   5492,
#   7130,
#   60664,
#   17432,
#   28396,
#   32086
#  ), genus = c(
#   16109,
#   47702,
#   168067,
#   7643,
#   72065,
#   122035,
#   156321,
#   69740
#  ), family = c(
#   16109,
#   47726,
#   183846,
#   72092,
#   122035,
#   156486,
#   69782
#  ))
# 
# I_04 <- list(
#  total = 1486238,
#  species = c(
#   1633,
#   44367,
#   15819,
#   21774,
#   107532,
#   25471,
#   50939,
#   78158
#  ), genus = c(
#   75121,
#   95955,
#   487478,
#   23492,
#   128593,
#   176800,
#   294170,
#   169474
#  ), family = c(
#   75121,
#   96007,
#   536884,
#   128626,
#   176805,
#   294500,
#   169571
#  ))
# 
# II_03 <- list(
#  total = 722657,
#  species = c(
#   383,
#   13072,
#   6970,
#   6582,
#   48061,
#   16771,
#   36617,
#   29889
#  ), genus = c(
#   13440,
#   28813,
#   201484,
#   7230,
#   59176,
#   105967,
#   214896,
#   77802
#  ), family = c(
#   13440,
#   28843,
#   218014,
#   59188,
#   105967,
#   215088,
#   77897
#  ))

# ## align
# I_03 <- list(
#  total = 549714,
#  species = c(
#   7046,
#   45204,
#   22333,
#   70206,
#   67399,
#   98605,
#   70868,
#   49420
#  ), genus = c(
#   13583,
#   47770,
#   80541,
#   72992,
#   67632,
#   100850,
#   113316,
#   52924
#  ), family = c(
#   13583,
#   47770,
#   153615,
#   67632,
#   100850,
#   113316,
#   52924
#  ))
# 
# I_04 <- list(
#  total = 1211528,
#  species = c(
#   31510,
#   91487,
#   63952,
#   209469,
#   122552,
#   141916,
#   125620,
#   119005
#  ), genus = c(
#   63536,
#   98681,
#   226654,
#   218195,
#   123158,
#   147450,
#   205347,
#   127607
#  ), family = c(
#   63536,
#   98681,
#   445214,
#   123158,
#   147450,
#   205347,
#   127607
#  ))
# 
# II_03 <- list(
#  total = 466470,
#  species = c(
#   2080,
#   24197,
#   19040,
#   61553,
#   54152,
#   72445,
#   39185,
#   29603
#  ), genus = c(
#   9950,
#   37216,
#   86546,
#   64035,
#   55637,
#   89505,
#   84915,
#   33565
#  ), family = c(
#   9950,
#   37216,
#   152067,
#   55637,
#   89520,
#   84915,
#   33565
#  ))

## otu
I_03 <- list(
 total = 497891,
 species = c(
   3,
   11318,
   922,
   234,
   23509,
   9220,
   22897,
   17707
 ), genus = c(
   11558,
   28150,
   139484,
   407,
   40685,
   88439,
   136728,
   47760
 ), family = c(
   11558,
   28167,
   140263,
   40685,
   88439,
   136735,
   47760
 ))

I_04 <- list(
 total = 1042934,
 species = c(
   107,
   5195,
   4688,
   8703,
   56749,
   23903,
   13266,
   55774
 ), genus = c(
   50732,
   43378,
   339450,
   8769,
   86540,
   122300,
   237465,
   115294
 ), family = c(
   50732,
   43378,
   367622,
   86540,
   122300,
   237469,
   115294
 ))

II_03 <- list(
 total = 433436,
 species = c(
   20,
   2161,
   1132,
   795,
   17493,
   16593,
   14276,
   15857
 ), genus = c(
   4474,
   22841,
   125761,
   995,
   34517,
   66196,
   135642,
   38000
 ), family = c(
   4474,
   22841,
   128735,
   34517,
   66196,
   135645,
   38006
 ))


reference <- rbind(
data.frame(
value = c(4.2,9.9,10.1,10.4,14.1,15.5,17.4,18.4,0),
L1 = "species"),
data.frame(
value = c(4.2,9.9,10.1,10.4,14.1,15.5,17.4,18.4,0),
L1 = "genus"),
data.frame(
value = c(4.2,9.9,20.5,14.1,15.5,17.4,18.4,0),
L1 = "family")
)
reference$type <- "reference"
reference$names <- factor(c(
 "P. aeruginosa",
 "E. faecalis",
 "E. coli",
 "S. enterica",
 "L. monocytogenes",
 "S. aureus",
 "B. subtilis",
 "L. fermentum",
 "other",
 
 "Pseudomonas",
 "Enterococcus",
 "Escherichia-Shigella",
 "Salmonella",
 "Listeria",
 "Staphylococcus",
 "Bacillus",
 "Lactobacillus",
 "other",
 
 "Pseudomonadaceae",
 "Enterococcaceae",
 "Enterobacteriaceae",
 "Listeriaceae",
 "Staphylococcaceae",
 "Bacillaceae",
 "Lactobacillaceae",
 "other"
), levels = c("P. aeruginosa","E. faecalis","E. coli","S. enterica","L. monocytogenes","S. aureus","B. subtilis","L. fermentum","Pseudomonas","Enterococcus","Escherichia-Shigella","Salmonella","Listeria","Staphylococcus","Bacillus","Lactobacillus","Pseudomonadaceae","Enterococcaceae","Enterobacteriaceae","Listeriaceae","Staphylococcaceae","Bacillaceae","Lactobacillaceae","other"
))
reference$relative <- 1

I_03[["species"]] <- (I_03[["species"]] / I_03[["total"]]) * 100
I_03[["species"]] <- c( I_03[["species"]], 100 - sum(I_03[["species"]]) )
I_03[["genus"]] <- (I_03[["genus"]] / I_03[["total"]]) * 100
I_03[["genus"]] <- c( I_03[["genus"]], 100 - sum(I_03[["genus"]]) )
I_03[["family"]] <- (I_03[["family"]] / I_03[["total"]]) * 100
I_03[["family"]] <- c( I_03[["family"]], 100 - sum(I_03[["family"]]) )
I_03 <- melt(I_03[2:4])
I_03$names <- reference$names
I_03$type = "run I ZyDNA"
I_03$relative <- I_03$value / reference$value


I_04[["species"]] <- (I_04[["species"]] / I_04[["total"]]) * 100
I_04[["species"]] <- c( I_04[["species"]], 100 - sum(I_04[["species"]]) )
I_04[["genus"]] <- (I_04[["genus"]] / I_04[["total"]]) * 100
I_04[["genus"]] <- c( I_04[["genus"]], 100 - sum(I_04[["genus"]]) )
I_04[["family"]] <- (I_04[["family"]] / I_04[["total"]]) * 100
I_04[["family"]] <- c( I_04[["family"]], 100 - sum(I_04[["family"]]) )
I_04 <- melt(I_04[2:4])
I_04$names <- reference$names
I_04$type = "run I ZyCell"
I_04$relative <- I_04$value / reference$value

II_03[["species"]] <- (II_03[["species"]] / II_03[["total"]]) * 100
II_03[["species"]] <- c( II_03[["species"]], 100 - sum(II_03[["species"]]) )
II_03[["genus"]] <- (II_03[["genus"]] / II_03[["total"]]) * 100
II_03[["genus"]] <- c( II_03[["genus"]], 100 - sum(II_03[["genus"]]) )
II_03[["family"]] <- (II_03[["family"]] / II_03[["total"]]) * 100
II_03[["family"]] <- c( II_03[["family"]], 100 - sum(II_03[["family"]]) )
II_03 <- melt(II_03[2:4])
II_03$names <- reference$names
II_03$type = "run II ZyDNA"
II_03$relative <- II_03$value / reference$value

df_all <- rbind(reference,I_03,I_04,II_03)
df_all$type <- factor(df_all$type, levels = c("reference","run I ZyDNA","run I ZyCell","run II ZyDNA"))

df_all

p_S <- ggplot(df_all[df_all$L1=="species",]) + theme_bw() + 
 geom_bar( aes( x = type, y = value, fill = names ), stat = "identity", color = "black" ) + 
 scale_fill_manual( values = c( brewer.pal(8, "Dark2"), "white" ) ) + 
 theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
 ylab("relative abundance") + labs(fill="species")

ggsave(filename = "species.png", plot = p_S, width = 10, height = 9, units = "cm")

p_G <- ggplot(df_all[df_all$L1=="genus",]) + theme_bw() + 
 geom_bar( aes( x = type, y = value, fill = names ), stat = "identity", color = "black" ) + 
 scale_fill_manual( values = c( brewer.pal(8, "Dark2"), "white" ) ) + 
 theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
 ylab("relative abundance") + labs(fill="genus")

ggsave(filename = "genus.png", plot = p_G, width = 10, height = 9, units = "cm")

p_F <- ggplot(df_all[df_all$L1=="family",]) + theme_bw() + 
 geom_bar( aes( x = type, y = value, fill = names ), stat = "identity", color = "black" ) + 
 scale_fill_manual( values = c( brewer.pal(8, "Dark2")[c(1:3,5:8)], "white" ) ) + 
 theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
 ylab("relative abundance") + labs(fill="family")

ggsave(filename = "family.png", plot = p_F, width = 10, height = 9, units = "cm")


df_all$relative <- log10(df_all$relative)

p_S <- ggplot(df_all[df_all$L1=="species" & df_all$names!="other" & df_all$type!="reference" ,]) + theme_bw() +
 geom_point( aes( y = relative, x = names, color = names), size = 2 ) + 
 geom_segment( aes( y = relative, yend = 0, x = names, xend = names, color = names) ) + 
 scale_color_manual( values = brewer.pal(8, "Dark2") ) + 
 geom_text(aes(y=0,x="P. aeruginosa",label="")) +
 theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
 ylab("abundance deviation [log10]") + labs(color="species") +
 facet_grid( type ~ . )

ggsave(filename = "species_dev.png", plot = p_S, width = 10, height = 9, units = "cm")

p_G <- ggplot(df_all[df_all$L1=="genus" & df_all$names!="other" & df_all$type!="reference" ,]) + theme_bw() +
 geom_point( aes( y = relative, x = names, color = names), size = 2 ) + 
 geom_segment( aes( y = relative, yend = 0, x = names, xend = names, color = names) ) + 
 scale_color_manual( values = brewer.pal(8, "Dark2") ) + 
 geom_text(aes(y=0,x="Pseudomonas",label="")) +
 theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
 ylab("abundance deviation [log10]") + labs(color="genus") +
 facet_grid( type ~ . )

ggsave(filename = "genus_dev.png", plot = p_G, width = 10, height = 9, units = "cm")

p_F <- ggplot(df_all[df_all$L1=="family" & df_all$names!="other" & df_all$type!="reference" ,]) + theme_bw() +
 geom_point( aes( y = relative, x = names, color = names), size = 2 ) + 
 geom_segment( aes( y = relative, yend = 0, x = names, xend = names, color = names) ) + 
 scale_color_manual( values = brewer.pal(8, "Dark2")[c(1:3,5:8)] ) + 
 geom_text(aes(y=0,x="Pseudomonadaceae",label="")) +
 theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
 ylab("abundance deviation [log10]") + labs(color="family") +
 facet_grid( type ~ . )

ggsave(filename = "family_dev.png", plot = p_F, width = 10, height = 9, units = "cm")
