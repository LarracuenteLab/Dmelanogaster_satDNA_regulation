#This script takes count and col info file output from summarize_multi_count_for_DESeq2.py, and run DESeq2 with them.
#author: Xiaolu Wei (xiaolu_wei@urmc.rochester.edu)

library("DESeq2")

args <- commandArgs(trailingOnly=TRUE)
countdata_file <- args[1]
coldata_file <- args[2]
out_file <- args[3]

#read in count data
countdata <- read.table(countdata_file, header = T, row.names = 1)
#filter out low count
countdata_keep <- countdata[rowSums(countdata)>10,]

#read in column info
coldata <- read.table(coldata_file, header = T, row.names = 1)
coldata$condition <- factor(coldata$condition, levels=c("wildtype","mutant"))

mydata <- DESeqDataSetFromMatrix(countData = countdata_keep,
                              colData = coldata,
                              design = ~ condition)

#test data from more than one paper
if (length(levels(coldata$study)) > 1) {
  design(mydata) <- formula(~ study + condition)
  print("Use both $study and $condition as covariant.")
} else{
  print("Use only $condition as variant.")
}

#run DESeq2
mydata <- DESeq(mydata)

result <- results(mydata, alpha=0.05)
mcols(result)$description

resOrdered <- result[order(result$pvalue),]

write.csv(as.data.frame(resOrdered), file=out_file)


