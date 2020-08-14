#piRNA_RPFRandom_BarplotZscore_BH.R
#This script takes the count (wc -l) file and draw barplot based on col1 and col2
#This is an adapted script for lendis of piTargets mapping, and calculate 41th Z score
#Version: Yu Sun, 2016-12-02

args <- commandArgs(TRUE)
Name <- args[1]

Data <- read.table(Name)

#Z-score
Rest <- c(Data$V2[1:9],Data$V2[11:length(Data$V2)])
Z_sore <- (Data$V2[10]-mean(Rest))/sd(Rest)
Z_sore

#Draw Lendis
File_pdf <- paste(Name,".pdf",sep="")
pdf(File_pdf,width=16,height=9)
barplot(
  Data$V2,names.arg=Data$V1,
  xlab="Transcipt",
  ylab='Count',
  main=Name,
  col="darkblue"
)
Subtitle <- paste("10th position Z score = ",Z_sore, sep = "")
mtext(Subtitle)
dev.off()
