### CD to directory with methylkit output images
setwd(
### Make design matrix
sample.id <- c("D1","D2","D3","D4","D5","D6","D8","D9","D10", "D11","D12","D13","D14","D15","D16","D17","D18","D19","D20")
treatment_n <- c(1,1,1,0,1,0,1,1,0,1,0,1,0,0,1,0,0,1,0)
treatment <- gsub(1,"sick",treatment_n) %>% gsub("0","healthy",.)
design.matrix <- bind_cols(ID=sample.id, treatment=treatment)

### Make PCA 
pca_5X <- percMethylation(meth.filtered_5X) %>% t(.) %>% prcomp(.,center=T, scale. = F)
percent_variance_pca_5X<-summary(pca_5X)$importance["Proportion of Variance",] * 100

pca_10X <- percMethylation(meth.filtered_10X) %>% t(.) %>% prcomp(.,center=T, scale. = F)
percent_variance_pca_10X<-summary(pca_10X)$importance["Proportion of Variance",] * 100

### plot PCA 
pca_plot_5X<-as_tibble(pca_5X$x) %>% bind_cols(design.matrix)
pca_plot_10X<-as_tibble(pca_10X$x) %>% bind_cols(design.matrix)

ggplot(pca_plot_5X,aes(x=PC1, y=PC2, color=treatment,label= ID)) +
  geom_point()  +
  xlab(label = paste("PC1 (",round(percent_variance_pca[1],1),"percent of the total variance)")) +
  ylab(label = paste("PC2 (",round(percent_variance_pca[2],1),"percent of the total variance)"))+
  ggtitle("Crappyfish PCA 5X")+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_blank())
       
