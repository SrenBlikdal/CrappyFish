library(methylKit)
library(tidyverse)
library(qqman)

#methylkitworkspace<-("/home/sren/Documents/Master_thesis/Crappyfish/methylkit_19samples_5X_210416.RData")
#methylkitworkspace<-("/home/sren/Documents/Master_thesis/Crappyfish/methylkit_19samples_10X_MN_210129.RData")
load(methylkitworkspace)
rm(myDiff_MN)
rm(meth.filtered)
#grep the chromosomes, to do, assign contigs to chr
myDiff_t<-as_tibble(myDiff_10X)
#rm(myDiff_MN)
#rm(meth.filtered)
#Find pvalue threshold for qvalue>0.01
#notsig<-myDiff_t %>% subset(qvalue>=0.01)
#myDiff_t<-myDiff_t %>% subset(qvalue<=0.01)

threshold_for_sig<--log10(min(notsig$pvalue))

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
  
rm(myDiff_t)

axisdf = SMP %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

ggplot(SMP, aes(x=BPcum, y=(SMvalue))) +
  
  # Show all points
  geom_point( aes(color=as.factor(chr)), alpha=0.6, size=1.0) +
  scale_color_manual(values = rep(c("orange", "black"), 600 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
  geom_hline(yintercept=threshold_for_sig, linetype="dashed", color = "red")+
  geom_hline(yintercept=-threshold_for_sig, linetype="dashed", color = "red")+
  xlab("Position in genome") +
  ylab("log10(p-value) for hypo and -log10(p-value) for hyper") +

  # Custom the theme:
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


