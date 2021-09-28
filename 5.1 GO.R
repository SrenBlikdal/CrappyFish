library(AnnotationHub)
library(AnnotationDbi)
library(readr)
library(dplyr)
library(clusterProfiler)

ah <- AnnotationHub()
# list OrgDb's available for Salmo salar
subset(ah, ah$species=="Salmo salar" & ah$rdataclass=="OrgDb")
# retrieve the Ssal OrgDb
SsalOrg <- ah[["AH80527"]]

DMC_X_PG<-read.csv("~/Documents/Master_thesis/Crappyfish/5X_DMC_X_PG.txt", header = F)
cov_X_PG<-read.csv("~/Documents/Master_thesis/Crappyfish/5X_cov_X_PG.txt", header = F)
dim(unique(DMC_X_PG))
dim(unique(cov_X_PG))
GO_PG_BF <- enrichGO(gene = as.character(unique(DMC_X_PG$V1)), 
                  pAdjustMethod = "fdr",
                  minGSSize = 1,
                  maxGSSize = 10000,
                  ont = "BP", # can only test one ontology at a time. BP = "Biological process"
                  universe = as.character(unique(cov_X_PG$V1)),
                  OrgDb = SsalOrg)

res_GO_PG_BF<-dplyr::select(GO_PG_BF@result,-geneID, -qvalue) %>% dplyr::filter(.,p.adjust<0.01)
DT::datatable((res_GO_PG_BF),rownames = F)
emapplot(GO_PG_BF, showCategory = nrow(res_GO_PG_BF))

GO_PG_MF <- enrichGO(gene = as.character(unique(DMC_X_PG$V1)), 
                     pAdjustMethod = "fdr",
                     minGSSize = 1,
                     maxGSSize = 10000,
                     ont = "MF", # can only test one ontology at a time. MF = "Molecular function"
                     universe = as.character(unique(cov_X_PG$V1)),
                     OrgDb = SsalOrg)

res_GO_PG_MF<-dplyr::select(GO_PG_MF@result,-geneID, -qvalue) %>% dplyr::filter(.,p.adjust<0.01)
DT::datatable((res_GO_PG_MF),rownames = F)
emapplot(GO_PG_MF,showCategory = nrow(res_GO_PG_MF))

GO_PG_CC <- enrichGO(gene = as.character(unique(DMC_X_PG$V1)), 
                     pAdjustMethod = "fdr",
                     minGSSize = 1,
                     maxGSSize = 10000,
                     ont = "CC", # can only test one ontology at a time. CC = "Cellular component"
                     universe = as.character(unique(cov_X_PG$V1)),
                     OrgDb = SsalOrg)

res_GO_PG_CC<-dplyr::select(GO_PG_CC@result,-geneID, -qvalue) %>% dplyr::filter(.,p.adjust<0.01)
DT::datatable((res_GO_PG_CC),rownames = F)
emapplot(GO_PG_CC, showCategory = nrow(res_GO_PG_CC))

resKEGG <- enrichKEGG(gene = DMC_X_PG$V1,
                      universe = as.character(cov_X_PG$V1),
                      pAdjustMethod = "BH",
                      organism = "sasa")

DT::datatable(dplyr::select(resKEGG@result,-geneID),rownames = F)
dim(unique(DMC_X_P))
dim(unique(cov_X_P))
DMC_X_P<-read.csv("~/Documents/Master_thesis/Crappyfish/5X_DMC_X_P.txt", header = F)
cov_X_P<-read.csv("~/Documents/Master_thesis/Crappyfish/5X_cov_X_P.txt", header = F)
GO_P_BP<- enrichGO(gene = as.character(unique(DMC_X_P$V1)), 
                  pAdjustMethod = "fdr",
                  minGSSize = 1,
                  maxGSSize = 10000,
                  ont = "BP", # can only test one ontology at a time. BP = "Biological process"
                  universe = as.character(unique(cov_X_P$V1)),
                  OrgDb = SsalOrg)

res_GO_P_BP<-dplyr::select(GO_P_BP@result,-geneID, -qvalue) %>% dplyr::filter(.,p.adjust<0.01)
DT::datatable((res_GO_P_BP),rownames = F)
emapplot(GO_P_BP)

GO_P_CC<- enrichGO(gene = as.character(unique(DMC_X_P$V1)), 
                   pAdjustMethod = "fdr",
                   minGSSize = 1,
                   maxGSSize = 10000,
                   ont = "CC", # can only test one ontology at a time. BP = "Biological process"
                   universe = as.character(unique(cov_X_P$V1)),
                   OrgDb = SsalOrg)

res_GO_P_CC<-dplyr::select(GO_P_CC@result,-geneID, -qvalue) %>% dplyr::filter(.,p.adjust<0.01)
DT::datatable((res_GO_P_CC),rownames = F)
emapplot(GO_P_CC)

GO_P_MF<- enrichGO(gene = as.character(unique(DMC_X_P$V1)), 
                   pAdjustMethod = "fdr",
                   minGSSize = 1,
                   maxGSSize = 10000,
                   ont = "MF", # can only test one ontology at a time. BP = "Biological process"
                   universe = as.character(unique(cov_X_P$V1)),
                   OrgDb = SsalOrg)

res_GO_P_MF<-dplyr::select(GO_P_MF@result,-geneID, -qvalue) %>% dplyr::filter(.,p.adjust<0.01)
DT::datatable((res_GO_P_MF),rownames = F)
emapplot(GO_P_MF)

dim(unique(cov_X_P))
DMR_X_P<-read.csv("~/Documents/Master_thesis/Crappyfish/5X_DMR_X_P.txt", header = F)
cov_X_P<-read.csv("~/Documents/Master_thesis/Crappyfish/5X_cov_X_P.txt", header = F)
GO_P_DMR_BP <- enrichGO(gene = as.character(unique(DMR_X_P$V1)), 
                 pAdjustMethod = "fdr",
                 minGSSize = 1,
                 maxGSSize = 10000,
                 ont = "BP", # can only test one ontology at a time. BP = "Biological process"
                 universe = as.character(unique(cov_X_P$V1)),
                 OrgDb = SsalOrg)

res_GO_P_DMR_BP<-dplyr::select(GO_P_DMR_BP@result,-geneID, -qvalue) %>% dplyr::filter(.,p.adjust<0.01)
DT::datatable((res_GO_P_DMR_BP),rownames = F)
emapplot(GO_P_DMR_BP)

GO_P_DMR_CC <- enrichGO(gene = as.character(unique(DMR_X_P$V1)), 
                        pAdjustMethod = "fdr",
                        minGSSize = 1,
                        maxGSSize = 10000,
                        ont = "CC", # can only test one ontology at a time. BP = "Biological process"
                        universe = as.character(unique(cov_X_P$V1)),
                        OrgDb = SsalOrg)

res_GO_P_DMR_CC<-dplyr::select(GO_P_DMR_CC@result,-geneID, -qvalue) %>% dplyr::filter(.,p.adjust<0.01)
DT::datatable((res_GO_P_DMR_CC),rownames = F)
emapplot(GO_P_DMR_CC)

GO_P_DMR_MF <- enrichGO(gene = as.character(unique(DMR_X_P$V1)), 
                        pAdjustMethod = "fdr",
                        minGSSize = 1,
                        maxGSSize = 10000,
                        ont = "MF", # can only test one ontology at a time. BP = "Biological process"
                        universe = as.character(unique(cov_X_P$V1)),
                        OrgDb = SsalOrg)

res_GO_P_DMR_MF<-dplyr::select(GO_P_DMR_MF@result,-geneID, -qvalue) %>% dplyr::filter(.,p.adjust<0.01)
DT::datatable((res_GO_P_DMR_MF),rownames = F)
emapplot(GO_P_DMR_MF)
