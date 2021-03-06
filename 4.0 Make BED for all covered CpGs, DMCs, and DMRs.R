library(tidyverse)
library(methylKit)

load("MethylKit_19samples_10X_5X_210901.RData")

### All sites covered to 5X
myDiff_5X_t <- as_tibble(myDiff_5X)
write.table(myDiff_5X_t, file = "/home/projects/ku-cbd/data/HoloFish/morten_fish/Epigenome/Novogene_CrappyFish_20samples_24Sep20_Analyses/5X_cov_CpGs.tsv", row.names=FALSE, sep="\t")

### All 5X DMCs 
DMCs_5X_t <- as_tibble(myDiff_5X) %>% filter(.,qvalue<0.01)
write.table(DMCs_5X_t, file = "/home/projects/ku-cbd/data/HoloFish/morten_fish/Epigenome/Novogene_CrappyFish_20samples_24Sep20_Analyses/5X_DMCs.tsv", row.names=FALSE, sep="\t")

### All 5X DMRs 
DMRs_5X_t <- as_tibble(mysigdmr_5X)
write.table(DMRs_5X_t, file = "/home/projects/ku-cbd/data/HoloFish/morten_fish/Epigenome/Novogene_CrappyFish_20samples_24Sep20_Analyses/5X_DMRs.tsv", row.names=FALSE, sep="\t")

### Convert to bedfiles in terminal 
sed 's/"//g' 5X_cov_CpGs.tsv | awk -v OFS='\t' '{print $1, $2, $3}' > 5X_cov_CpGs.bed
sed 's/"//g' 5X_DMCs.tsv | awk -v OFS='\t' '{print $1, $2, $3}' > 5X_DMCs.bed
sed 's/"//g' 5X_DMRs.tsv | awk -v OFS='\t' '{print $1, $2, $3}' > 5X_DMRs.bed
##Remove header with nano 
