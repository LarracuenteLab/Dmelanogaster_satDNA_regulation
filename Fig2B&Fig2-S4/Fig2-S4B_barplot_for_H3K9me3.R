#author: Xiaolu Wei (xiaolu_wei@urmc.rochester.edu)

library(ggplot2)
library(dplyr)

count_all=read.table("summary_H3K9me3.txt",header=FALSE)
summary(count_all)

count <-count_all[count_all$V2!='353bp_SAT' & count_all$V2!='356bp_SAT',]

count$V2 <- factor(count$V2, 
                   levels=c("Rsp_SAT","260bp_SAT","359bp_SAT","piRNA_cluster_20A","piRNA_cluster_flamenco","piRNA_cluster_42AB","piRNA_cluster_80F","piRNA_cluster_38C1","piRNA_cluster_38C2"), 
                   labels = c("Rsp","260bp","359bp","20A","flamenco","42AB","80F","38C1","38C2"))


#get statistics summary
count_sum <- ddply(count, c("V2"), summarise,
                   N    = length(V3),
                   mean = mean(V3),
                   sd   = sd(V3),
                   se   = sd / sqrt(N) )

pdf("barplot_ovary_H3K9me3_chip.pdf",width=5,height=2.5)
ggplot(count_sum, aes(V2, mean))+
  geom_bar(stat="identity", position = "dodge", width=0.8, fill="skyblue3", color="black") + 
  #error bars
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, size=0.3, position=position_dodge(.9))+
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


