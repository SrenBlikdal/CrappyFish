library(tidyverse)
#Make tibble with number of overlaps between all CpG sites, all CpG sites covered in the analysis, DMC, hyper DMCs and hypo DMCs and the genetic features . 

type<-rep(c('1Kpromotor', 'gene','250promotor', '6KPromotor'),5)
context<-c(rep("All_CpG",4),rep("All_CpG_cov",4),rep("DMCs",4),rep("DMC_hyper",4),rep("DMC_hypo",4))
#imported values from 4_gene_annotation_terminal
cov<-c(420871, 11928830, 230762, 2272035, 363882, 9082103, 156039, 2144700,772,9837,189,3186,507,6528,121,2062,265,3309,68,1124)

context<-c(rep("All_CpG",4),rep("All_cov_CpG",4),rep("DMCs",4),rep("4hypo DMRs",4))

dat=tibble(type,context,cov)
ggplot(dat,aes(x=context,y=cov, fill=type)) +
    geom_bar(position="fill", stat="identity")
chisq.test(matrix(c(363882,16615963,772,19013),ncol = 2))
chisq.test(matrix(c(363882,16615963,772,19013),ncol = 2))
