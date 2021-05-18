#author: Xiaolu Wei (xiaolu_wei@urmc.rochester.edu)

library(ggplot2)
library(plyr)

count_all=read.table("summary_H3K9me3_wt.txt",header=FALSE)
summary(count_all)

#boxplot
count <- count_all

count$V2 <- factor(count$V2, 
                   levels=c("Rsp_SAT","1pt688_SAT","1pt688.2L_2_260bp_locus","1pt688.3L_locus","1pt688.Contig101_locus","1pt688.Contig9_locus","piRNA_cluster_20A","piRNA_cluster_flamenco","piRNA_cluster_42AB","piRNA_cluster_80F","piRNA_cluster_38C1","piRNA_cluster_38C2","euchromatin"), 
                   )

for ( i in 1:nrow(count)) {
  if(count[i,2]=="1pt688.2L_2_260bp_locus" || count[i,2]=="1pt688.3L_locus" || count[i,2]=="1pt688.Contig101_locus" || count[i,2]=="1pt688.Contig9_locus" ) {
    count[i,2]="1pt688_SAT"
  }
}

count$V2 <- factor(count$V2, 
                   levels=c("Rsp_SAT","1pt688_SAT","piRNA_cluster_20A","piRNA_cluster_flamenco","piRNA_cluster_42AB","piRNA_cluster_80F","piRNA_cluster_38C1","piRNA_cluster_38C2","euchromatin"), 
                   labels = c("Rsp","1pt688","20A","flamenco","42AB","80F","38C1","38C2","euchromatin"))

pdf("boxplot_ovary_H3K9me3_chip_1pt688_locus.pdf",width=7,height=3)
ggplot(count, aes(V2, V3))+
  geom_boxplot(show.legend=T, color="black", fill="skyblue3") + 
  xlab("")+ 
  ylab("H3K9me3 ChIP/Input")+
  theme(panel.background = element_rect(colour = "black", fill = "white"))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank())+
  geom_hline(yintercept=1, size=0.25,colour="dark grey")+
  theme(axis.text.y = element_text(colour = "black", size=10))+
  #face="italic"
  theme(axis.text.x = element_text(colour = "black", size=10))+
  theme(axis.ticks = element_line(colour = "black"))+
  theme(axis.line = element_line(size=0.5))+
  theme(axis.title = element_text(size=10))+
  theme(legend.title = element_blank())
dev.off()

