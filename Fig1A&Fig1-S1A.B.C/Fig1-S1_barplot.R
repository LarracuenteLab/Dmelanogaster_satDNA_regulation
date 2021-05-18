#@author: Xiaolu Wei (xiaolu_wei@urmc.rochester.edu)

library(ggplot2)
library(dplyr)

#dev.stages

count_all=read.csv("BDGP_dev.stage_sat.csv",header=FALSE)
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

pdf("totalRNA_dev.stage.pdf",width=5,height=2)
ggplot(count_sum, aes(factor(V2,levels=c("0-2hr","12-14hr","22-24hr","L2","WPP","WPP_24hr","WPP_3d","male_1d","fem_1d","male_30d","fem_30d")), mean))+
  geom_bar(stat="identity", position = "dodge", width=0.7, fill="#56B4E9", color="black") + 
  #error bars
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, size=0.3, position=position_dodge(.9))+
  xlab("")+ 
  ylab("RPM")+
  theme(panel.background = element_rect(colour = "black", fill = "white"))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank())+
  geom_hline(yintercept=0, size=0.25,colour="grey")+
  theme(axis.text.y = element_text(colour = "black", size=6))+
  theme(axis.text.x = element_text(colour = "black", size=6))+
  theme(axis.ticks = element_line(colour = "black"))+
  theme(axis.line = element_line(size=0.5))+
  theme(axis.title = element_text(size=6))+
  theme(legend.title = element_blank())
dev.off()


##################

#tissues

count_all=read.csv("BDGP_tissue_total_sat.csv", header=FALSE)
summary(count_all)
#For 1.688
count=count_all[count_all$V1=='1.688',]
#OR
#For Rsp
count=count_all[count_all$V1=='Rsp',]

#get statistics summary
count_sum <- ddply(count, c("V2"), summarise,
                   N    = length(V3),
                   mean = mean(V3),
                   sd   = sd(V3),
                   se   = sd / sqrt(N) )

pdf("totalRNA_tissue.pdf",width=5,height=3)
ggplot(count_sum, aes(factor(V2,levels=c("ovary","testes","head","carcass")), mean))+
  geom_bar(stat="identity", position = "dodge", width=0.7, fill="light green", color="black") + 
  #error bars
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, size=0.3, position=position_dodge(.9))+
  xlab("")+ 
  ylab("RPM")+
  theme(panel.background = element_rect(colour = "black", fill = "white"))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank())+
  geom_hline(yintercept=0, size=0.25,colour="grey")+
  theme(axis.text.y = element_text(colour = "black", size=12))+
  theme(axis.text.x = element_text(colour = "black", size=12))+
  theme(axis.ticks = element_line(colour = "black"))+
  theme(axis.line = element_line(size=0.5))+
  theme(axis.title = element_text(size=12))+
  theme(legend.title = element_blank())
dev.off()

