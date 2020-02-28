# 
# rule deseq2:
#     input:
#         expand("{tmp}01_processed_data/03_kmer_mapping/{run}/{barc}/{reference}_{reftype}/Species.counttaxlist", 
#         tmp = config["experiments"]["tmp"], run = RUNS, barc = SAMPLES.keys(), reference = config["reference"]["source"], reftype = config["reference"]["type"])
#     output:
#         "doesntexist"
#     shell:
#         "echo {input}"
# 
# 
# 













# > ls()
# [1] "count_tab"       "sample_info_tab" "tax_tab"
# 
# > str(count_tab)
# 'data.frame':   2498 obs. of  16 variables:
#  $ BW1  : int  0 0 0 0 0 0 0 0 0 0 ...
#  $ BW2  : int  0 0 1 0 0 0 0 1 0 0 ...
#  $ R10  : int  0 459 44 392 0 76 160 11 11 5 ...
#  $ R11BF: int  1762 1 0 0 0 4 0 3 3 0 ...
#  $ R11  : int  4 130 26 236 0 82 48 0 8 4 ...
#  $ R12  : int  4 78 148 6 0 74 117 77 76 71 ...
#  $ R1A  : int  0 43 204 2 0 64 71 115 183 70 ...
#  $ R1B  : int  0 38 153 0 0 74 81 80 258 102 ...
#  $ R2   : int  3 90 154 38 0 115 125 104 72 118 ...
#  $ R3   : int  0 148 68 64 0 128 90 94 41 116 ...
#  $ R4   : int  5 159 113 62 0 155 166 213 111 284 ...
#  $ R5   : int  7 71 140 56 0 108 78 107 150 154 ...
#  $ R6   : int  4 8 193 0 0 122 29 77 129 87 ...
#  $ R7   : int  2 0 5 0 0 81 0 0 49 48 ...
#  $ R8   : int  2 178 39 260 0 38 123 98 0 0 ...
#  $ R9   : int  0 159 19 160 1199 56 29 134 11 5 ...
#  > head(count_tab)
#        BW1 BW2 R10 R11BF R11 R12 R1A R1B  R2  R3  R4  R5  R6 R7  R8   R9
#  ASV_1   0   0   0  1762   4   4   0   0   3   0   5   7   4  2   2    0
#  ASV_2   0   0 459     1 130  78  43  38  90 148 159  71   8  0 178  159
#  ASV_3   0   1  44     0  26 148 204 153 154  68 113 140 193  5  39   19
#  ASV_4   0   0 392     0 236   6   2   0  38  64  62  56   0  0 260  160
#  ASV_5   0   0   0     0   0   0   0   0   0   0   0   0   0  0   0 1199
#  ASV_6   0   0  76     4  82  74  64  74 115 128 155 108 122 81  38   56
# 
# > str(tax_tab)
#  chr [1:2498, 1:6] "Bacteria" "Bacteria" "Bacteria" "Bacteria" "Bacteria" ...
#  - attr(*, "dimnames")=List of 2
#   ..$ : chr [1:2498] "ASV_1" "ASV_2" "ASV_3" "ASV_4" ...
#   ..$ : chr [1:6] "Kingdom" "Phylum" "Class" "Order" ...
#  > head(tax_tab)
#        Kingdom    Phylum           Class                 Order
#  ASV_1 "Bacteria" "Proteobacteria" "Alphaproteobacteria" "Kordiimonadales"
#  ASV_2 "Bacteria" "Nitrospirae"    "Nitrospira"          "Nitrospirales"
#  ASV_3 "Bacteria" "Proteobacteria" "Alphaproteobacteria" "Rhodovibrionales"
#  ASV_4 "Bacteria" "Proteobacteria" "Gammaproteobacteria" "pItb-vmat-80"
#  ASV_5 "Bacteria" "Proteobacteria" "Gammaproteobacteria" "UBA10353_marine_group"
#  ASV_6 "Bacteria" "Proteobacteria" "Gammaproteobacteria" "Steroidobacterales"
#        Family           Genus
#  ASV_1 NA               NA
#  ASV_2 "Nitrospiraceae" "Nitrospira"
#  ASV_3 "Kiloniellaceae" NA
#  ASV_4 NA               NA
#  ASV_5 NA               NA
#  ASV_6 "Woeseiaceae"    "Woeseia"
# 
# > str(sample_info_tab)
# 'data.frame':   16 obs. of  4 variables:
#  $ temp : num  2 2 13.7 7.3 7.3 NA 8.6 8.6 8.6 12.7 ...
#  $ type : Factor w/ 3 levels "biofilm","rock",..: 3 3 2 1 2 2 2 2 2 2 ...
#  $ char : Factor w/ 5 levels "altered","biofilm",..: 5 5 4 2 4 1 1 1 1 1 ...
#  $ color: chr  "blue" "blue" "black" "darkgreen" ...
