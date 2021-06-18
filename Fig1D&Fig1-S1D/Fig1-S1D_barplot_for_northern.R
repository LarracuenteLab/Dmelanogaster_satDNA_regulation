############
#R script to make Figure 1-supplement 1D. Plot showing the Northern signal 
#quantification of Rsp transcript abundance  in different fly strains in 
#Figure 1C.
#@author: Xiaolu Wei (xiaolu_wei@urmc.rochester.edu)
###########

library(ggplot2)

count=read.csv(file.choose(),header=FALSE)
summary(count)

count=read.table(file.choose(),header=FALSE)
summary(count)

pdf("Rsp.northern.pdf",width=5,height=3)
ggplot(count, aes(factor(V1, levels=c("ZW144","Ral357","Iso1","Ral380","It pk cnbw")),V2)) +
  geom_bar(stat="identity", position="identity", width=0.5,show.legend=F, fill="darkorchid4") + 
  xlab("")+
  ylab("Relative Rsp hybridization")+
  theme(panel.background = element_rect(colour = "black", fill = "white"))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank())+
  geom_hline(yintercept=0, size=0.25,colour="grey")+
  theme(axis.text.y = element_text(colour = "black", size=10))+
  theme(axis.text.x = element_text(colour = "black", size=10))+
  theme(axis.ticks = element_line(colour = "black"))+
  theme(axis.line = element_line(size=0.5))+
  theme(axis.title = element_text(size=10))+
  theme(legend.title = element_blank() )
dev.off()

# correlation
a <- c(238,640,1100,2339,4086)
b <- c(0.966552, 6.798309, 7.934042, 8.456499, 17.28229)
cor(a, b, method="pearson")
cor.test(a, b, alternative = "two.sided", method = "pearson")

