#author: Xiaolu Wei (xiaolu_wei@urmc.rochester.edu)

library(ggplot2)
library(dplyr)

count_all=read.table("summary_ovary_Akkouche_H3K9me3_piwiemKD_fold_chagne.txt",header=F)
summary(count_all)

count <- count_all
count$V2 <- factor(count$V2, 
                   levels=c("Rsp_SAT","1pt688_SAT","piRNA_cluster_20A","piRNA_cluster_flamenco","piRNA_cluster_42AB","piRNA_cluster_80F","piRNA_cluster_38C1","piRNA_cluster_38C2"), 
                   labels = c("Rsp","1.688","20A","flamenco","42AB","80F","38C1","38C2"))


#boxplot
pdf("boxplot_ovary_Akkouche_H3K9me3_piwiemKD_fold_change.pdf",width=5, height=2)
ggplot(count, aes(V2, V3))+
  geom_boxplot(show.legend=F, color="black", fill="orange") + 
  xlab("")+ 
  ylab("Log2FC of enrichment")+
  theme(panel.background = element_rect(colour = "white", fill = "white"))+
  theme(panel.grid.major = element_blank())+
  geom_hline(yintercept=0, size=0.25,colour="grey")+
  theme(axis.text.y = element_text(colour = "black", size=10))+
  theme(axis.text.x = element_text(colour = "black", size=10, face="italic"))+
  theme(axis.ticks = element_line(colour = "black"))+
  theme(axis.line = element_line(size=0.5))+
  theme(axis.title = element_text(size=10))+
  theme(legend.title = element_blank())
dev.off()


#t test
t.test(count[count$V2=="Rsp",]$V3, mu=0)
t.test(count[count$V2=="1.688",]$V3, mu=0)

t.test(count[count$V2=="20A",]$V3, mu=0)
t.test(count[count$V2=="flamenco",]$V3, mu=0)

t.test(count[count$V2=="42AB",]$V3, mu=0)
t.test(count[count$V2=="80F",]$V3, mu=0)
t.test(count[count$V2=="38C1",]$V3, mu=0)
t.test(count[count$V2=="38C2",]$V3, mu=0)

p.adjust(c(0.001839,0.0005507,
           0.07337,0.006191,
           0.006005,3.585e-05,0.02291,0.05989), method="BH")

