#author: Xiaolu Wei (xiaolu_wei@urmc.rochester.edu)

library(ggplot2)
library(dplyr)

#depth

count_all=read.table("SRR1187952.260bp.depth",header=FALSE)
summary(count_all)
count <- count_all

#Rsp (exclude TEs)
count <- count[count$V2<12722 | count$V2>16829, ]
count <- count[count$V2<18608 | count$V2>22480, ]
count <- count[count$V2<23844 | count$V2>27254, ]
count <- count[count$V2<111760 | count$V2>113086, ]
count <- count[count$V2<129632 | count$V2>133739, ]
count <- count[count$V2<135518 | count$V2>135900, ]
count <- count[count$V2<137158 | count$V2>140605, ]
count <- count[count$V2<141969 | count$V2>145841, ]
count <- count[count$V2<147618 | count$V2>151674, ]
count <- count[count$V2<154081 | count$V2>158673, ]
count <- count[count$V2<159347 | count$V2>165491, ]

#OR
#260bp (exclude TEs)
count <- count[count$V2<431869 | count$V2>437016, ]
count <- count[count$V2<442185 | count$V2>447332, ]
count <- count[count$V2<459709 | count$V2>460136, ]

pdf("ovary_260bp.pdf",width=5,height=2)
ggplot(count, aes(x=factor(V2), y=V3, color=V4))+
  geom_bar(stat="identity", position = "identity", width=0.3)+
  xlab("position (kb)")+ 
  ylab("read depth")+
  #change based on data
#  scale_y_continuous(breaks=c(-10,-5,0,5,10),limits=c(-10,10))+
  #Rsp
#  scale_x_discrete(breaks=c("50000","100000","150000"), labels=c("50","100","150"))+
  #OR
  #260bp, limits=c(402701,460225)
  scale_x_discrete(breaks=c("410000","440000"), labels=c("410","440"))+
  scale_color_manual(values=c("orange","blue"), guide=FALSE)+
  theme(panel.background = element_rect(colour = "white", fill = "white"))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank())+
  geom_hline(yintercept=0, size=0.5,colour="black")+
  theme(axis.text = element_text(colour = "black", size=5))+
  theme(axis.ticks = element_line(colour = "black"))+
  theme(axis.line = element_line(size=0.5))+
  theme(axis.title = element_text(size=8))
dev.off()

##############################

###count###

count_all=read.csv("stranded_count.csv",header=TRUE)
summary(count_all)

count <- count_all[count_all$type=="minus_all" | count_all$type=="plus_all", ]
#OR
count <- count_all[count_all$type=="minus_unique" | count_all$type=="plus_unique", ]

#get statistics summary
count_sum <- ddply(count, c("repeat.","type"), summarise,
                   N    = length(ratio),
                   mean = mean(ratio),
                   sd   = sd(ratio),
                   se   = sd / sqrt(N) )
summary(count_sum)

pdf("stranded_count.pdf",width=3.5,height=2)
ggplot(count_sum, aes(x=factor(type,levels=c("plus_all","minus_all","plus_unique","minus_unique")), y=mean, fill=factor(repeat.,levels=c("Rsp","1.688"))))+
  geom_bar(stat="identity", position = "dodge", width=0.8)+
  xlab("")+ 
  ylab("% of reads")+
  scale_fill_manual(values=c("red2","gold2"))+
  #error bars
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.3, size=0.3, position=position_dodge(.8))+
  theme(panel.background = element_rect(colour = "white", fill = "white"))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank())+
  geom_hline(yintercept=0, size=0.25,colour="dark grey")+
  theme(axis.text.y = element_text(colour = "black", size=10))+
  #face="italic"
  theme(axis.text.x = element_text(colour = "black", size=10))+
  theme(axis.ticks = element_line(colour = "black"))+
  theme(axis.line = element_line(size=0.5))+
  theme(axis.title = element_text(size=10))+
  theme(legend.title = element_blank())
dev.off()


