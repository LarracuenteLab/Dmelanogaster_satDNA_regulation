##########
#R script to make Figure 2B and Figure 2-figure supplement 5A. Figure 2B showing 
#Rhino ChIP-seq enrichment scores for satDNAs, uni-strand (uni) piRNA clusters, 
#dual-strand (dual) and euchromatin (eu). P-values are estimated by pairwise 
#t-tests with FDR correction (Benjamini 1995). Figure 2-figure supplement 5A showing
#the enrichment scores for each satDNA and piRNA cluster.
#author: Xiaolu Wei (xiaolu_wei@urmc.rochester.edu)
########

library(ggplot2)
library(plyr)

count_all=read.table("summary_ovary_rhino.txt",header=FALSE)
summary(count_all)

#####
#Figure 2-figure supplement 5A
#####
count <- count_all
count$V2 <- factor(count$V2, 
                   levels=c("Rsp_SAT","1pt688_SAT","1pt688.2L_2_260bp_locus","1pt688.3L_locus","1pt688.Contig101_locus","1pt688.Contig9_locus","piRNA_cluster_20A","piRNA_cluster_flamenco","piRNA_cluster_42AB","piRNA_cluster_80F","piRNA_cluster_38C1","piRNA_cluster_38C2","euchromatin"), 
)

for ( i in 1:nrow(count)) {
if(count[i,2]=="1pt688.2L_2_260bp_locus" || count[i,2]=="1pt688.3L_locus" || count[i,2]=="1pt688.Contig101_locus" || count[i,2]=="1pt688.Contig9_locus" ) {
  count[i,2]="1pt688_SAT"
}
}

count$V2 <- factor(count$V2, 
                   levels=c("Rsp_SAT","1pt688_SAT","piRNA_cluster_20A","piRNA_cluster_flamenco","piRNA_cluster_42AB","piRNA_cluster_80F","piRNA_cluster_38C1","piRNA_cluster_38C2","euchromatin"), 
                   labels = c("Rsp","1pt688","20A","flamenco","42AB","80F","38C1","38C2","euchromatin"))

pdf("boxplot_ovary_rhino_chip_1pt688_locus.pdf",width=7,height=3)
ggplot(count, aes(V2, V3))+
  geom_boxplot(show.legend=T, color="black", fill="brown") + 
  xlab("")+ 
  ylab("Rhino ChIP/Input")+
  theme(panel.background = element_rect(colour = "black", fill = "white"))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank())+
  geom_hline(yintercept=1, size=0.25,colour="dark grey")+
  theme(axis.text.y = element_text(colour = "black", size=10))+
  #face="italic"
  theme(axis.text.x = element_text(colour = "black", size=10))+
  theme(axis.ticks = element_line(colour = "black"))+
  theme(axis.line = element_line(size=0.5))+
  theme(axis.title = element_text(size=10))+
  theme(legend.title = element_blank())
dev.off()


#####
#Figure 2B 
#####
##combined plot
count <- count_all
count$V4 <- ifelse(count$V2=="piRNA_cluster_flamenco" | count$V2=="piRNA_cluster_20A", "uni", 
                   ifelse(count$V2=="piRNA_cluster_42AB" | count$V2=="piRNA_cluster_80F" |count$V2=="piRNA_cluster_38C1" | count$V2=="piRNA_cluster_38C2", 
                           "dual", ifelse(count$V2=="euchromatin", "eu", "SAT" )))
count$V4 <- factor(count$V4, 
                  levels=c("SAT", "uni","dual","eu"))

#get statistics summary including standard error
count_sum <- ddply(count, c("V4"), summarise,
                   N    = length(V3),
                   mean = mean(V3),
                   sd   = sd(V3),
                   se   = sd / sqrt(N) )

pdf("barplot_ovary_rhino_chip_combined_1pt688_locus_uniq.pdf",width=2.5,height=2)
ggplot(count_sum, aes(V4, mean))+
  geom_bar(stat="identity", position = "dodge", width=0.8, fill="brown", color="black") + 
  #error bars
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, size=0.3, position=position_dodge(.9))+
  #color
  #scale_fill_manual(values=c("#56B4E9"))+
  xlab("")+ 
  ylab("Rhino ChIP/Input")+
  theme(panel.background = element_rect(colour = "white", fill = "white"))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank())+
  geom_hline(yintercept=1, size=0.25,colour="dark grey")+
  theme(axis.text.y = element_text(colour = "black", size=10))+
  #face="italic"
  theme(axis.text.x = element_text(colour = "black", size=10))+
  theme(axis.ticks = element_line(colour = "black"))+
  theme(axis.line = element_line(size=0.5))+
  theme(axis.title = element_text(size=10))+
  theme(legend.title = element_blank())
dev.off()


#significant test

#Bartlett test of homogeneity of varianc
bartlett.test(V3 ~ V4, data=count)  #p-value = 2.021e-12
##pairwise t.test
pairwise.t.test(x=count$V3, g=count$V4, p.adjust.method="BH", pool.sd=FALSE)

