### CD to directory with methylkit output images
setwd()
load()
### Make design matrix
sample.id <- c("D1","D2","D3","D4","D5","D6","D8","D9","D10", "D11","D12","D13","D14","D15","D16","D17","D18","D19","D20")
treatment_n <- c(1,1,1,0,1,0,1,1,0,1,0,1,0,0,1,0,0,1,0)
treatment <- gsub(1,"Sick",treatment_n) %>% gsub("0","Healthy",.)
design.matrix <- bind_cols(ID=sample.id, treatment=treatment)

### Make PCA 
pca_5X <- percMethylation(meth.filtered_5X) %>% t(.) %>% prcomp(.,center=T, scale. = F)
percent_variance_pca_5X<-summary(pca_5X)$importance["Proportion of Variance",] * 100

pca_10X <- percMethylation(meth.filtered_10X) %>% t(.) %>% prcomp(.,center=T, scale. = F)
percent_variance_pca_10X<-summary(pca_10X)$importance["Proportion of Variance",] * 100

### plot PCA 
pca_plot_5X<-as_tibble(pca_5X$x) %>% bind_cols(design.matrix)
pca_plot_10X<-as_tibble(pca_10X$x) %>% bind_cols(design.matrix)

png(file="/home/projects/ku-cbd/data/HoloFish/morten_fish/Epigenome/Novogene_CrappyFish_20samples_24Sep20_Analyses/PCA_5X.png")
ggplot(pca_plot_5X,aes(x=PC1, y=PC2, color=treatment,label= ID)) +
  geom_point(size=3)  +
  geom_text_repel(size=7) +
  scale_color_manual(values=c("blue", "red"))+
  xlab(label = paste("PC1 (",round(percent_variance_pca_5X[1],1),"percent of the total variance)")) +
  ylab(label = paste("PC2 (",round(percent_variance_pca_5X[2],1),"percent of the total variance)"))+
  ggtitle("PCA \n 5X threshold")+
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")+
  theme_bw(base_size=16)+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(color=guide_legend("Disease status"))
  dev.off()
  
png(file="/home/projects/ku-cbd/data/HoloFish/morten_fish/Epigenome/Novogene_CrappyFish_20samples_24Sep20_Analyses/PCA_10X.png")
ggplot(pca_plot_10X,aes(x=PC1, y=PC2, color=treatment,label= ID)) +
  geom_point(size=3)  +
  geom_text_repel(size=7) +
  scale_color_manual(values=c("blue", "red"))+
  xlab(label = paste("PC1 (",round(percent_variance_pca_10X[1],1),"percent of the total variance)")) +
  ylab(label = paste("PC2 (",round(percent_variance_pca_10X[2],1),"percent of the total variance)"))+
  ggtitle("PCA \n 10X threshold")+
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")+
  theme_bw(base_size=16)+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(color=guide_legend("Disease status"))
  dev.off()
  
       
