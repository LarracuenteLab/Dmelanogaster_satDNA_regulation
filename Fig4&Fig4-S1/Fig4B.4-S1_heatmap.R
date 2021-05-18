#author: Xiaolu Wei (xiaolu_wei@urmc.rochester.edu)

library(ggplot2)
library(dplyr)

#total RNA
count_all=read.csv("summary-totalRNA-embryonic-piwiKD.csv",header=T)
count <- count_all
count$repeatID <- factor(count$repeatID, 
                         levels=c("Rsp_SAT","1pt688_SAT","piRNA_cluster_20A","piRNA_cluster_flamenco","piRNA_cluster_42AB","piRNA_cluster_80F","piRNA_cluster_38C1","piRNA_cluster_38C2"), 
                         labels = c("Rsp","1.688","20A","flamenco","42AB","80F","38C1","38C2"))
#heatmap
pdf("heatmap_ovary_Akkouche_totalRNA_piwiemKD_fold_change.pdf", width=7, height=1)
ggplot(count, aes(repeatID,protein)) + 
  geom_tile(aes(fill = log2FoldChange)) + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-2.6,2.6) ) +  
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) 
dev.off()


#small RNA
count_all=read.table("smallRNA_norm_miRNA_embryonic_piwiKD.log2FC.txt",header=F)
count <- count_all
count$V3 <- factor(count$V3, 
                         levels=c("Rsp_SAT","1pt688_SAT","piRNA_cluster_20A","piRNA_cluster_flamenco","piRNA_cluster_42AB","piRNA_cluster_80F","piRNA_cluster_38C1","piRNA_cluster_38C2"), 
                         labels = c("Rsp","1.688","20A","flamenco","42AB","80F","38C1","38C2"))
#heatmap
pdf("heatmap_ovary_Akkouche_smallRNA_piwiemKD_fold_change.pdf", width=7, height=2)
ggplot(count, aes(V3,V2)) + 
  geom_tile(aes(fill = V4)) + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-2.6,2.6) ) +  
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) 
dev.off()


