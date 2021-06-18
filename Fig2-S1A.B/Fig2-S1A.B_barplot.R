############
#R script to make Figure 2-figure supplement 1A-B. Plot showing small RNA levels 
#in RPM (reads per million) in different tissues for 1.688 and Rsp.
#author: Xiaolu Wei (xiaolu_wei@urmc.rochester.edu)
############

library(ggplot2)
library(dplyr)

count_all=read.csv("smallRNA_tissues.csv",header=FALSE)
summary(count_all)
#for 1.688
count=count_all[count_all$V1=='1.688',]
#OR
#for Rsp
count=count_all[count_all$V1=='Rsp',]

#get statistics summary
count_sum <- ddply(count, c("V2"), summarise,
                   N    = length(V3),
                   mean = mean(V3),
                   sd   = sd(V3),
                   se   = sd / sqrt(N) )

pdf("smallRNA_tissue.pdf",width=5,height=3)
ggplot(count_sum, aes(factor(V2,levels=c("ovary","testis","head","whole body")), mean))+
  geom_bar(stat="identity", position = "dodge", width=0.7, fill="light green", color="black") + 
  #error bars
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, size=0.3, position=position_dodge(.9))+
  xlab("")+ 
  ylab("RPM")+
  theme(panel.background = element_rect(colour = "black", fill = "white"))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank())+
  geom_hline(yintercept=0, size=0.25,colour="grey")+
  theme(axis.text.y = element_text(colour = "black", size=15))+
  theme(axis.text.x = element_text(colour = "black", size=15))+
  theme(axis.ticks = element_line(colour = "black"))+
  theme(axis.line = element_line(size=0.5))+
  theme(axis.title = element_text(size=15))+
  theme(legend.title = element_blank())
dev.off()


