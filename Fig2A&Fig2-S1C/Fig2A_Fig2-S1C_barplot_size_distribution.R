#author: Xiaolu Wei (xiaolu_wei@urmc.rochester.edu)

library(ggplot2)
library(dplyr)

#with error bars

count=read.table("mix_sex_Rsp.size", header=FALSE)
count=read.table("mix_sex_1.688.size", header=FALSE)

summary(count)

#get statistics summary including standard error
count_sum <- ddply(count, c("V1"), summarise,
                   N    = length(V2),
                   mean = mean(V2),
                   sd   = sd(V2),
                   se   = sd / sqrt(N) )

pdf("mix_sex_size_distribute_Rsp_errorbar.pdf",width=10,height=5)
ggplot(count_sum, aes(factor(V1), mean))+
  geom_bar(stat="identity", position = "dodge", fill='#56B4E9', color="black") + 
  #error bars
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.3, size=0.5, position=position_dodge(.9))+
  xlab("Small RNA size (nt)")+ 
  ylab("% of total mapped reads")+
  theme(panel.background = element_rect(colour = "black", fill = "white"))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank())+
  geom_hline(yintercept=0, size=0.25,colour="grey")+
  theme(axis.text.y = element_text(colour = "black", size=15))+
  theme(axis.text.x = element_text(colour = "black", size=15))+
  theme(axis.ticks = element_line(colour = "black"))+
  theme(axis.line = element_line(size=0.5))+
  theme(axis.title = element_text(size=18))+
  theme(legend.title = element_blank())
dev.off()


