# R script to find differentially expressed genes by R package DESeq2
library("DESeq2")
library("limma")
library("dplyr")

d<-read.table("filtered_samples.30bin.B_N.remove_outlier.counts",header = T,sep = "\t")
counts <- d[,2:42]
rownames(counts)<-d[,1]
keep <- rowSums(counts>=1)>=9
counts <- counts[keep,]
group<-c(rep("HCC",30),rep("Normal",11))

colData <- data.frame(row.names = colnames(counts),group_list=group)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData,design = ~group_list)
dds2 <- DESeq(dds)


norm_value <- rlogTransformation(dds2)
exprSet=assay(norm_value)
write.csv(exprSet,file = "DESeq2_norm_value.csv",row.names = T)

resultsNames(dds2)
results <- results(dds2,contrast = c("group_list","HCC","Normal"))
write.table(results,"DESeq2_results.csv",row.names = T,sep = ",")

uniq=na.omit(results)
uniq$up_down=ifelse(uniq$log2FoldChange>0,"up","down")
write.csv(uniq,file="DESeq2_results_uniq.csv",row.names = T)

padj_0.05=uniq[uniq$padj<0.05,]
write.csv(padj_0.05,file = "DESeq2_results_padj_0.05.csv")

up_gene<-padj_0.05[padj_0.05$log2FoldChange>=1,]
down_gene<-padj_0.05[padj_0.05$log2FoldChange<=-1,]
dim(up_gene)
dim(down_gene)
write.csv(up_gene,"DESeq2_up.csv")
write.csv(down_gene,"DESeq2_down.csv")
