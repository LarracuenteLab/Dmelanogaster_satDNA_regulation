##########
#R script to make Figure 3-figure supplement 2. Plot showing the quantification 
#of total RNA abundance for piRNA clusters in RDC (A) and Rhi/Moon (B) mutant 
#ovaries relative to wild type (log2 fold change of 1-kb windows).
#author: Xiaolu Wei (xiaolu_wei@urmc.rochester.edu)
########

library(ggplot2)

#read in "count_unique_1kbwindow_new.csv.summary_RDC_piRNAclusters" for (A) and
#"count_unique_1kbwindow.csv.summary_Rhi_Moon_piRNAclusters" for (B).
count=read.csv("count_unique_1kbwindow_new.csv.summary_RDC_piRNAclusters", header=T)

count$newrepeatID <- factor(count$newrepeatID, 
                            levels=c("piRNA_cluster_20A","piRNA_cluster_flamenco","piRNA_cluster_42AB","piRNA_cluster_80F","piRNA_cluster_38C1","piRNA_cluster_38C2"), 
                            labels = c("20A","flamenco","42AB","80F","38C1","38C2"))

# for (A) RDC
count$protein <- factor(count$protein,
                        levels = c("Rhino_1", "Deadlock","Cutoff"),
                        labels = c("Rhino", "Deadlock","Cutoff"))

#OR
#for (B) Rhi_Moon
count$protein <- factor(count$protein,
                        levels = c("Rhino_2", "Moonshiner"),
                        labels = c("Rhino", "Moonshiner"))

#boxplot
pdf("boxplot_totalRNA_1kbwindow_unique_RDC.pdf",width=7.5, height=2)
ggplot(count, aes(newrepeatID, log2FoldChange, fill=protein))+
  geom_boxplot(show.legend=T, color="black") + 
  #color 
  #"gold","light blue"
  #"gold","yellow3","yellowgreen"
  scale_fill_manual(values=c("gold","light blue"))+
  xlab("")+ 
  ylab("Log2FC")+
  theme(panel.background = element_rect(colour = "white", fill = "white"))+
  theme(panel.grid.major = element_blank())+
  #theme(panel.grid.minor = element_blank())+
  geom_hline(yintercept=0, size=0.25,colour="grey")+
  theme(axis.text.y = element_text(colour = "black", size=10))+
  theme(axis.text.x = element_text(colour = "black", size=10, face="italic"))+
  theme(axis.ticks = element_line(colour = "black"))+
  theme(axis.line = element_line(size=0.5))+
  theme(axis.title = element_text(size=10))+
  theme(legend.title = element_blank())
dev.off()


