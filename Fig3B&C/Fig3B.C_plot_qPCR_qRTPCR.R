##########
#R script to make Figure 3B and C. Figure 3B showing qPCR estimate of Rsp copy 
#number in wild types and mutants. Figure 3C showing qRT-PCR estimate of Rsp 
#transcript level in mutants compared to wild types. P-values are calculated by
#Student’s t-test.
#author: Xiaolu Wei (xiaolu_wei@urmc.rochester.edu)
########

library(ggplot2)
library(dplyr)


#Figure 3B (bar plot)

count=read.csv("moon.rhi.mutant.DNA.qPCR.csv",header=FALSE)
summary(count)

#get statistics summary including standard error
count_sum <- ddply(count, c("V1"), summarise,
                   N    = length(V2),
                   mean = mean(V2),
                   sd   = sd(V2),
                   se   = sd / sqrt(N) )

pdf("Rsp.moon.rhi.mutants.DNA_1.20.2020.pdf",width=4, height=2)
ggplot(count_sum, aes(factor(V1,levels=c("w1118","w1-OregonR","moon mutant","rhi mutant"), labels=c("w1118","w1/OreR","moon mutant","rhi mutant")), mean ))+
  geom_bar(stat="identity", position = "identity", width=0.8,show.legend=F, color="black", fill="light green") + 
  #error bars
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, size=0.5, position=position_dodge(.9))+
  xlab("")+ 
  ylab("Rsp copy number")+
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


#Figure 3C (boxplot)

count=read.csv("moon.rhi.mutant.RNA.qRT-PCR.csv",header=FALSE)
summary(count)

#get statistics summary including standard error
count_sum <- ddply(count, c("V1"), summarise,
                   N    = length(V2),
                   mean = mean(V2),
                   sd   = sd(V2),
                   se   = sd / sqrt(N) )

pdf("Rsp.moon.rhi.mutants.expression.normalized.by.DNA.copy.number_1.20.2020.pdf",width=3, height=2)
ggplot(count, aes(factor(V1,levels=c("moon mutant","rhi mutant")), V2 ))+
  geom_boxplot(show.legend=F, color="black", fill="light blue") + 
  xlab("")+ 
  ylab("delta delta Ct")+
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

#paired t test
t.test(count[count$V1=="moon mutant",]$V2, mu=0)
t.test(count[count$V1=="rhi mutant",]$V2, mu=0)


