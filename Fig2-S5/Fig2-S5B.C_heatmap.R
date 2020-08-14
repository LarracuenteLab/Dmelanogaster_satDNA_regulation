#author: Xiaolu Wei (xiaolu_wei@urmc.rochester.edu)

library(ggplot2)
library(dplyr)

#total RNA
count_all=read.csv("summary-totalRNA-embryonic-piwiKD_combined.csv",header=T)

#small RNA
count_all=read.table("smallRNA_norm_miRNA_embryonic_piwiKD.log2FC_combined_3.txt",header=F)

summary(count_all)

count <-count_all[count_all$repeatID!='353bp_SAT' & count_all$repeatID!='356bp_SAT',]

count$repeatID <- factor(count$repeatID, 
                         levels=c("Rsp_SAT","260bp_SAT","359bp_SAT","piRNA_cluster_20A","piRNA_cluster_flamenco","piRNA_cluster_42AB","piRNA_cluster_80F","piRNA_cluster_38C1","piRNA_cluster_38C2"), 
                         labels = c("Rsp","260bp","359bp","20A","flamenco","42AB","80F","38C1","38C2"))

#heatmap
pdf("heatmap_ovary_Akkouche_totalRNA_piwiemKD_fold_change.pdf", width=7, height=1)
ggplot(count, aes(repeatID,protein)) + 
  geom_tile(aes(fill = log2FoldChange)) + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-2.6,2.6) ) +  
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) 
dev.off()
