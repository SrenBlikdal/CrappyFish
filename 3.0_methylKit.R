#!/usr/bin/Rscript

###
library(BiocManager)
library(memoise)
library(sessioninfo)
library(devtools)
library(methylKit)

setwd("/home/projects/ku-cbd/data/HoloFish/morten_fish/Epigenome/Novogene_CrappyFish_20samples_24Sep20_Analyses/5.1_bismark")

file.list<-list(
"D1.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
"D2.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
"D3.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
"D4.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
"D5.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
"D6.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
"D8.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
"D9.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
"D10.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
"D11.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
"D12.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
"D13.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
"D14.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
"D15.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
"D16.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
"D17.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
"D18.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
"D19.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
"D20.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz") 

my_obj<-methRead(file.list, sample.id = list("D1","D2","D3","D4","D5","D6","D8","D9","D10", "D11","D12","D13","D14","D15","D16","D17","D18","D19","D20"), assembly = "ICSASG_v2", pipeline = "bismarkCoverage", treatment= c(1,1,1,0,1,0,1,1,0,1,0,1,0,0,1,0,0,1,0), context= "CpG", mincov=5)
#treatment is meta data from Davides records:0 is healthy/Mycoplasma, 1 correspond to sick/Avivibro. 
#context =CpG, means we only look at CpG sites 

filtered.my_obj_1X=filterByCoverage(my_obj,lo.count=1,lo.perc=NULL, hi.count=NULL,hi.perc=99.9)
filtered.my_obj_5X=filterByCoverage(my_obj,lo.count=5,lo.perc=NULL, hi.count=NULL,hi.perc=99.9)
filtered.my_obj_10X=filterByCoverage(my_obj,lo.count=10,lo.perc=NULL, hi.count=NULL,hi.perc=99.9)

meth.filtered_1X<-unite(filtered.my_obj_5X, destrand=F)
meth.filtered_5X<-unite(filtered.my_obj_5X, destrand=F)
meth.filtered_10X<-unite(filtered.my_obj_10X, destrand=F)

print(dim(meth.filtered_5X))
print(dim(meth.filtered_10X))

myDiff_5X<-calculateDiffMeth(meth.filtered_5X, mc.cores=40)
myDiff_10X<-calculateDiffMeth(meth.filtered_10X, mc.cores=40)

mydmr_5X=edmr(myDiff_5X, mode=1, ACF=TRUE)
mysigdmr_5X=filter.dmr(mydmr, DMR.qvalue = 0.01, mean.meth.diff = 0, num.CpGs = 3, num.DMCs = 1)

mydmr_10X=edmr(myDiff_10X, mode=1, ACF=TRUE)
mysigdmr_10X=filter.dmr(mydmr, DMR.qvalue = 0.01, mean.meth.diff = 0, num.CpGs = 3, num.DMCs = 1)

#clean
rm(my_obj)
rm(filtered.my_obj_5X)
rm(filtered.my_obj_10X)

save.image(file ="MethylKit_19samples_210824.RData")
q()
