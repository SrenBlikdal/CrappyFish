library(tidyverse)
library(DescTools)
type<-rep(c("250promotor","1Kpromotor", "6Kpromotor","exon","intron","intergenic"),4)
cov<-c(404991,859448,4290210,2845060,16088490,23815697,150833,344881,1907524,1056458,6542839,6574559,181,738,2725,977,6786,7402,19,53,238,78,460,527)

context<-c(rep("All_CpG",6),rep("All_CpG_cov",6),rep("DMCs",6),rep("DMRs",6))

dat=tibble(type,context,cov)
pos<-dat %>% group_by(context) %>% summarise(sum_cov=sum(cov))
  
ggplot(dat,aes(x=context,y=cov,label=cov, fill=factor(type,levels = c("250promotor","1Kpromotor", "6Kpromotor","exon","intron","intergenic")))) +
  geom_col(position="fill")+
  geom_text(aes(label=stat(y),group =context), stat = "summary", fun =sum, nudge_y= (-(pos$sum_cov)+1.02)) +
  scale_fill_manual("Genetic feature", values = c("250promotor" = "blue4", "1Kpromotor" = "blue1", "6Kpromotor" = "deepskyblue1", "exon" ="green4", "intron"= "green1", "intergenic"= "indianred3"))+
  xlab("CpG Type")+
  ylab("Proportion of overlap")+
  theme_bw()

#G-test combined for all feature types testing:
#H0 = No difference in distribution of DMCs compared to the distribution of covered CpG-sites 
#H1 = Difference in distribution of DMCs compared to the distribution of covered CpG-sites
GTest(matrix(c(cov[7:12],cov[13:18]),nrow = 2, byrow = T))

#G-test for PP
GTest(matrix(c(cov[7],sum(cov[8:9],cov[10:13]),cov[13],sum(cov[14:18])),nrow = 2, byrow = T))
#G-test for MP
GTest(matrix(c(cov[8],sum(cov[7],cov[9:12]),cov[14],sum(cov[13],cov[15:18])),nrow = 2, byrow = T))
#G-test for DP
GTest(matrix(c(cov[9],sum(cov[7:8],cov[10:12]),cov[15],sum(cov[13:14],cov[16:18])),nrow = 2, byrow = T))
#G-test for exons
GTest(matrix(c(cov[10],sum(cov[7:9],cov[11:12]),cov[16],sum(cov[13:15],cov[17:18])),nrow = 2, byrow = T))
#G-test for introns
GTest(matrix(c(cov[11],sum(cov[7:10],cov[12]),cov[17],sum(cov[13:16],cov[18])),nrow = 2, byrow = T))
#G-test for intergenic regions
GTest(matrix(c(cov[12],sum(cov[7:11]),cov[18],sum(cov[13:17])),nrow = 2, byrow = T))
