library("tidyverse")

load("/home/projects/ku-cbd/data/HoloFish/morten_fish/Epigenome/Novogene_CrappyFish_20samples_24Sep20_Analyses/MethylKit_19samples_10X_5X_210901.RData")

DMCs_5X <- as_tibble(myDiff_5X) %>% filter(qvalue<0.01)
DMCs_10X <- as_tibble(myDiff_10X) %>% filter(qvalue<0.01)

DMCs_5X %>% filter(meth.diff>0) %>% dim(.)
DMCs_5X %>% filter(meth.diff<0) %>% dim(.)

DMCs_10X %>% filter(meth.diff>0) %>% dim(.)
DMCs_10X %>% filter(meth.diff<0) %>% dim(.)

DMRs_5X <- as_tibble(mysigdmr_5X)
DMRs_10X <- as_tibble(mysigdmr_10X)

DMRs_5X %>% filter(mean.meth.diff>0) %>% dim(.)
DMRs_5X %>% filter(mean.meth.diff<0) %>% dim(.)

DMRs_10X %>% filter(mean.meth.diff>0) %>% dim(.)
DMRs_10X %>% filter(mean.meth.diff<0) %>% dim(.)
