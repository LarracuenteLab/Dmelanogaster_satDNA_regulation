##########
#R script to make Figure 3 and Figure 3-figure supplement 1 and 3. Plot showing
#the quantification of small RNA abundance from piRNA pathway mutants for satDNAs 
#and piRNA clusters, normalized to the miRNA, or flamenco piRNA cluster.
#author: Xiaolu Wei (xiaolu_wei@urmc.rochester.edu)
########

library(ggplot2)

#read in one of the following files

#norm by miRNA
count<-read.table("summary-rhino-cuff-deadlock-moon.log2FC_norm_miRNA.txt",header=F)
count<-read.table("summary-primary-production.log2FC_norm_miRNA.txt",header=F)
count<-read.table("summary-pingpong.log2FC_norm_miRNA.txt",header=F)

count$V3 <- factor(count$V3, 
                   levels=c("Rsp_SAT","1pt688_SAT","piRNA_cluster_20A","piRNA_cluster_flamenco","piRNA_cluster_42AB","piRNA_cluster_80F","piRNA_cluster_38C1","piRNA_cluster_38C2"), 
                   labels = c("Rsp","1.688","20A","flamenco","42AB","80F","38C1","38C2"))

##OR
#norm by flam
count_all<-read.table("summary-rhino-cuff-deadlock-moon.log2FC_norm_flam.txt",header=F)
count_all<-read.table("summary-pingpong.log2FC_norm_flam.txt",header=F)

count <-count_all[count_all$V3!='piRNA_cluster_flamenco', ]
count$V3 <- factor(count$V3, 
                   levels=c("Rsp_SAT","1pt688_SAT","piRNA_cluster_20A","piRNA_cluster_flamenco","piRNA_cluster_42AB","piRNA_cluster_80F","piRNA_cluster_38C1","piRNA_cluster_38C2"), 
                   labels = c("Rsp","1.688","20A","flamenco","42AB","80F","38C1","38C2"))

#order the rows in the plots
########
#1 RDC
count$V2 <- factor(count$V2, 
                   levels=rev(c("Rhino_1","Rhino_2","Rhino_3","Rhino_4","Deadlock_1","Deadlock_2","Cutoff_1","Cutoff_2","Cutoff_3","Moonshiner")))

########
#2 primary
count$V2 <- factor(count$V2, 
                   levels=rev(c("Zucchini_1","Zucchini_2","Zucchini_3","Armitage_1","Armitage_2","Armitage_4","Gasz_1","Vreteno_1","Vreteno_2","Shutdown_1","Shutdown_2","Shutdown_3")))

########
#3 pingpong
count$V2 <- factor(count$V2, 
                   levels=rev(c("Ago3","Aub_1","Aub_2","Aub_3","Krimper_1","Krimper_2","SpnE_1","SpnE_2","Vasa_1","Vasa_2","UAP56","UAP56-oxidized","Piwi_1","Piwi_2","Piwi_3")))


#plot
pdf("heatmap_log2FC.pdf", width=10, height=5)
ggplot(count, aes(V3,V2)) + 
  geom_tile(aes(fill = V4)) + 
  xlab("") +
  ylab("") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits=c(-9.5,9.5), breaks=c(-9,-6,-3,0,3,6,9), name="log2 fold change") +  
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) 
dev.off()


