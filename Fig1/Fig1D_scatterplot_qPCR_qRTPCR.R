#@author: Xiaolu Wei (xiaolu_wei@urmc.rochester.edu)

library(ggplot2)
library(dplyr)

#diff.strain

#with error bars
count_DNA=read.csv("diff.strain.DNA.qPCR.csv",header=FALSE)
summary(count_DNA)
#get statistics summary including standard error
count_sum_DNA <- ddply(count_DNA, c("V1"), summarise,
                   N    = length(V2),
                   mean = mean(V2),
                   sd   = sd(V2),
                   se   = sd / sqrt(N) )

count_RNA=read.csv("diff.strain.RNA.qRT-PCR.csv",header=FALSE)
summary(count_RNA)
#get statistics summary including standard error
count_sum_RNA <- ddply(count_RNA, c("V1"), summarise,
                       N    = length(V2),
                       mean = mean(V2),
                       sd   = sd(V2),
                       se   = sd / sqrt(N) )

count_sum <- full_join(count_sum_DNA, count_sum_RNA, by="V1")
                    
pdf("diff.strain_scatter_plot.pdf",width=4,height=3, useDingbats=FALSE)
ggplot(count_sum, aes(mean.x, mean.y))+
  #point
  geom_point(color="dark blue", size=2) + 
  #error bars
  geom_errorbar(aes(ymin=mean.y-se.y, ymax=mean.y+se.y), width=0, size=0.2)+
  geom_errorbarh(aes(xmin = mean.x-se.x, xmax = mean.x+se.x), height=0, size=0.2)+
  geom_abline(intercept=4.340e-04, slope=2.942e-06, linetype="dotted", color="brown", size=0.5)+
  xlab("Rsp repeat number")+ 
  ylab("Rsp expression level")+
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

lm(mean.y ~ mean.x, data = count_sum)

