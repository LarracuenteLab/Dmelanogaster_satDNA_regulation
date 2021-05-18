#author: Xiaolu Wei (xiaolu_wei@urmc.rochester.edu)

library(ggplot2)
library(dplyr)

# change 1.688 to Rsp for plotting Rsp

count_ovary=read.table("ovary_1.688.size",header=FALSE)
count_ovary$tissue="ovary"

count_testes=read.table("testes_1.688.size",header=FALSE)
count_testes$tissue="testis"

count=rbind(count_ovary, count_testes)
count$tissue <- factor(count$tissue, 
                       levels = c("ovary", "testis"))
summary(count)

#get statistics summary including standard error
count_sum <- ddply(count, c("V1","tissue"), summarise,
                   N    = length(V2),
                   mean = mean(V2),
                   sd   = sd(V2),
                   se   = sd / sqrt(N) )

pdf("size_distribute_1pt688.pdf",width=6,height=2.5)
ggplot(count_sum, aes(factor(V1), mean, fill=tissue))+
  geom_bar(stat = "identity", position = "dodge", color="black") + 
  # , guide=FALSE
  #error bars
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.3, size=0.5, position=position_dodge(.9))+
  #color
  scale_fill_manual(values=c("lightblue", "royalblue4"))+
  xlab("Small RNA size (nt)")+ 
  ylab("% of total mapped reads")+
  theme(panel.background = element_rect(colour = "black", fill = "white"))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank())+
  geom_hline(yintercept=0, size=0.25,colour="grey")+
  theme(axis.text.y = element_text(colour = "black", size=10))+
  theme(axis.text.x = element_text(colour = "black", size=10))+
  theme(axis.ticks = element_line(colour = "black"))+
  theme(axis.line = element_line(size=0.5))+
  theme(axis.title = element_text(size=10))+
  theme(legend.title = element_blank())
dev.off()
