#author: Xiaolu Wei (xiaolu_wei@urmc.rochester.edu)

library(ggplot2)
library(ggseqlogo)
library(dplyr)

reads=read.table("ovary_1.688_24nt.txt",header=FALSE)
reads_vector <- as.vector(reads$V1)

pdf("ovary_1.688_24nt.pdf",width=4,height=2)
ggplot()+
  geom_logo(data=reads_vector, seq_type="rna")+
  theme(panel.background = element_rect(colour = "white", fill = "white"))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank())+
  geom_hline(yintercept=0, size=0.25,colour="grey")+
  theme(axis.text.y = element_text(colour = "black", size=6))+
  theme(axis.text.x = element_text(colour = "black", size=6))+
  theme(axis.ticks = element_line(colour = "black"))+
  theme(axis.line = element_line(size=0.3))+
  theme(axis.title = element_text(size=10))+
  theme(legend.title = element_blank())
dev.off()

