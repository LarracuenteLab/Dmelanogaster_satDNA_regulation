###########
#R script to make Figure 2–figure supplement 2 & 3. Plot showing read distribution
#along repeat consensus sequences.
#author: Xiaolu Wei (xiaolu_wei@urmc.rochester.edu)
##########

library(ggplot2)
library(dplyr)

###
#read in files for genomic varants
###

count=read.table("dmel_scaffold2_plus0310_Rsp_vs_monomer.distrib",header=TRUE)
repeatID="Rsp"
count$condition <- "genomic"
summary(count)

#Rsp 3L locus
count_2=read.table("dmel_scaffold2_plus0310_Rsp_3L_vs_monomer.distrib",header=TRUE)
repeatID="Rsp.with.Rsp_3L"
count_combined <- count
count_combined$condition <- "genomic"
count_combined$only3L <- count_2$total_count
count_combined$genome <- count$total_count


###
#read in files for Rhino ChIP-seq results
###

#Rsp as example
repeatID="Rsp"
suffix="perc.distrib"

treatment_rep1 <- read.table(paste("SRR5803102",repeatID,"vs_monomer",suffix, sep="_"), header=TRUE)
treatment_rep1$condition <- "treatment_rep1"

treatment_rep2 <- read.table(paste("SRR5806896",repeatID,"vs_monomer",suffix, sep="_"), header=TRUE)
treatment_rep2$condition <- "treatment_rep2"

treatment_rep3 <- read.table(paste("SRR1024325",repeatID,"vs_monomer",suffix, sep="_"), header=TRUE)
treatment_rep3$condition <- "treatment_rep3"

treatment_rep4 <- read.table(paste("SRR1024307",repeatID,"vs_monomer",suffix, sep="_"), header=TRUE)
treatment_rep4$condition <- "treatment_rep4"

#combine
count <- bind_rows(treatment_rep1, treatment_rep2, treatment_rep3, treatment_rep4)
#format condition
count$condition <- factor(count$condition, 
                          levels=c("treatment_rep1", "treatment_rep2", "treatment_rep3", "treatment_rep4"), 
                          labels=c("rhino_chip-1", "rhino_chip-2", "rhino_chip-3", "rhino_chip-4"))
summary(count)
count <- count[count$condition=="rhino_chip-1" | count$condition=="rhino_chip-2" | count$condition=="rhino_chip-3"| count$condition=="rhino_chip-4",]


###
#read in files for smallRNA-seq results
###

#Rsp as example
repeatID="Rsp"

rep1 <- read.table(paste("SRR5445299",repeatID,"vs_monomer.distrib", sep="_"), header=TRUE)
rep1$condition <- "rep1"

rep2 <- read.table(paste("SRR1187947",repeatID,"vs_monomer.distrib", sep="_"), header=TRUE)
rep2$condition <- "rep2"

rep3 <- read.table(paste("SRR5803096",repeatID,"vs_monomer.distrib", sep="_"), header=TRUE)
rep3$condition <- "rep3"

rep4 <- read.table(paste("SRR060645",repeatID,"vs_monomer.distrib", sep="_"), header=TRUE)
rep4$condition <- "rep4"

rep5 <- read.table(paste("SRR2046463",repeatID,"vs_monomer.distrib", sep="_"), header=TRUE)
rep5$condition <- "rep5"

rep6 <- read.table(paste("Iso1_R1",repeatID,"vs_monomer.distrib", sep="_"), header=TRUE)
rep6$condition <- "rep6"

#combine
count <- bind_rows(rep1, rep2, rep3, rep4, rep5, rep6)
#format condition
count$condition <- factor(count$condition, 
                          levels=c("rep1", "rep2", "rep3", "rep4", "rep5", "rep6"), 
                          labels=c("ovary-1", "ovary-2", "ovary-3", "testis-1", "testis-2", "testis-3"))

summary(count)

count <- count[count$condition=="testis-1" | count$condition=="testis-2" | count$condition=="testis-3",]
#OR
#count <- count[count$condition=="ovary-1" | count$condition=="ovary-2" | count$condition=="ovary-3",]


###
#PLOT HERE
###

outfile <- paste(repeatID, "lineplot.pdf", sep="_")
pdf(outfile, width=4, height=2)
ggplot(count, aes(position, total_count, color=condition))+
  geom_line() +
  #color
  #"black"
  #"burlywood1","burlywood2","burlywood3","burlywood4"
  #"orange", "violet", "red",
  scale_color_manual(values=c("green","#56B4E9","blue"))+
  xlab(repeatID)+ 
  ylab("%reads pile up")+
  theme(panel.background = element_rect(colour = "black", fill = "white"))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank())+
  geom_hline(yintercept=0, size=0.25,colour="grey")+
  theme(axis.text.y = element_text(colour = "black", size=8))+
  theme(axis.text.x = element_text(colour = "black", size=8))+
  theme(axis.ticks = element_line(colour = "black"))+
  theme(axis.line = element_line(size=0.5))+
  theme(axis.title = element_text(size=8))+
  theme(legend.title = element_blank())
dev.off()


#plot genomic Rsp distribution in solidline, and 3L Rsp distribution in dashed line
outfile <- paste(repeatID, "lineplot.pdf", sep="_")
pdf(outfile, width=4, height=2)
ggplot(count_combined, aes(x=position))+
  geom_line(aes(y=genome, color=condition)) +
  geom_line(aes(y=only3L, color=condition), linetype=2) +
  #color
  scale_color_manual(values=c("black"))+
  xlab(repeatID)+ 
  ylab("%reads pile up")+
  theme(panel.background = element_rect(colour = "black", fill = "white"))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank())+
  geom_hline(yintercept=0, size=0.25,colour="grey")+
  theme(axis.text.y = element_text(colour = "black", size=8))+
  theme(axis.text.x = element_text(colour = "black", size=8))+
  theme(axis.ticks = element_line(colour = "black"))+
  theme(axis.line = element_line(size=0.5))+
  theme(axis.title = element_text(size=8))+
  theme(legend.title = element_blank())
dev.off()


