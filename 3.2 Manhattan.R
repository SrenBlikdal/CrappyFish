library(methylKit)
library(tidyverse)
library(qqman)
load("/home/projects/ku-cbd/data/HoloFish/morten_fish/Epigenome/Novogene_CrappyFish_20samples_24Sep20_Analyses/MethylKit_19samples_10X_5X_210901.RData")

###5X sunset manhattan
myDiff_t <- as_tibble(myDiff_5X)

#Find pvalue threshold for qvalue>0.01

notsig <- myDiff_t %>% subset(qvalue>=0.01)

threshold_for_sig <- log10(min(notsig$pvalue))

SMP <- myDiff_t %>%
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(end)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(myDiff_t, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, end) %>%
  mutate( BPcum=end+tot) %>%
  mutate(hypo= ifelse((meth.diff<0),-1,1)) %>%
  mutate(SMvalue=hypo*(-log10(pvalue))) %>%
  dplyr::mutate_if(.,is.character,stringr::str_replace_all, pattern = "NW_*.*", replacement = "Scaffolds") %>%
  mutate(BPcum= if_else(chr == "Scaffolds", BPcum+50000000, BPcum))


axisdf = SMP %>% group_by(chr) %>% summarize(center=(max(BPcum) + min(BPcum) ) / 2 )
png(file="/home/projects/ku-cbd/data/HoloFish/morten_fish/Epigenome/Novogene_CrappyFish_20samples_24Sep20_Analyses/SNP_5X.png", width = 1000, height = 480)
ggplot(SMP, aes(x=BPcum, y=(SMvalue))) +
  
  # Show all points
  geom_point(aes(color=as.factor(chr)), alpha=0.6, size=1.0) +
  scale_color_manual(values = rep(c("orange", "black"),600)) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
  geom_hline(yintercept=threshold_for_sig, linetype="dashed", color = "red")+
  geom_hline(yintercept=-threshold_for_sig, linetype="dashed", color = "red")+
  xlab("Position in genome") +
  ylab("-log10(p-value)") +
  ylim(-30,30)+
  ggtitle("Sunset Manhattan Plot \n 5X threshold ")+

  # Custom the theme:
  theme_bw(base_size=16) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(hjust = 0.5))
dev.off()
###10X sunset manhattan
myDiff_t <- as_tibble(myDiff_10X)

#Find pvalue threshold for qvalue>0.01

notsig <- myDiff_t %>% subset(qvalue>=0.01)

threshold_for_sig <- log10(min(notsig$pvalue))

SMP <- myDiff_t %>%
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(end)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(myDiff_t, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, end) %>%
  mutate( BPcum=end+tot) %>%
  mutate(hypo= ifelse((meth.diff<0),-1,1)) %>%
  mutate(SMvalue=hypo*(-log10(pvalue))) %>%
  dplyr::mutate_if(.,is.character,stringr::str_replace_all, pattern = "NW_*.*", replacement = "Scaffolds") %>%
  mutate(BPcum= if_else(chr == "Scaffolds", BPcum+50000000, BPcum))


axisdf = SMP %>% group_by(chr) %>% summarize(center=(max(BPcum) + min(BPcum) ) / 2 )
png(file="/home/projects/ku-cbd/data/HoloFish/morten_fish/Epigenome/Novogene_CrappyFish_20samples_24Sep20_Analyses/SMP_10X.png", width = 1000, height = 480)
ggplot(SMP, aes(x=BPcum, y=(SMvalue))) +
  
  # Show all points
  geom_point(aes(color=as.factor(chr)), alpha=0.6, size=1.0) +
  scale_color_manual(values = rep(c("orange", "black"),600)) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
  geom_hline(yintercept=threshold_for_sig, linetype="dashed", color = "red")+
  geom_hline(yintercept=-threshold_for_sig, linetype="dashed", color = "red")+
  xlab("Position in genome") +
  ylab("-log10(p-value)") +
  ylim(-30,30)+
  ggtitle("Sunset Manhattan Plot \n 10X threshold ")+

  # Custom the theme:
  theme_bw(base_size=16) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(hjust = 0.5))
dev.off()

png(file="/home/projects/ku-cbd/data/HoloFish/morten_fish/Epigenome/Novogene_CrappyFish_20samples_24Sep20_Analyses/QQ_5X.png")
qq(myDiff_5X$pvalue, main = "Q-Q plot of EWAS p-values  \n 5X ", xlim = c(0, 30), ylim = c(0, 30))
dev.off()

png(file="/home/projects/ku-cbd/data/HoloFish/morten_fish/Epigenome/Novogene_CrappyFish_20samples_24Sep20_Analyses/QQ_10X.png")
qq(myDiff_10X$pvalue, main = "Q-Q plot of EWAS p-values  \n 10X ", xlim = c(0, 30), ylim = c(0, 30))
dev.off()
