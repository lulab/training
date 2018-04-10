# R script to find differentially expressed genes by R package edgeR.

library(limma)
library(edgeR)

d<-read.table("filtered_samples.30bin.B_N.remove_outlier.counts",header = T, sep = "\t")
counts <- d[,2:42]
rownames(counts)<-d[,1]
keep <- rowSums(counts>=1)>=9
counts <- counts[keep,]
group<-c(rep("HCC",30),rep("Normal",11))

DE_list <- DGEList(counts = counts,group = group)
DE_list <- calcNormFactors(DE_list)

design <- model.matrix(~group)
design <- model.matrix(~0+group,data = DE_list$samples)
colnames(design) <- levels(DE_list$samples$group)

DE_list <- estimateDisp(DE_list,design)


fit <- glmFit(DE_list,design)
lrt <- glmLRT(fit,contrast = c(1,-1))
FDR <- p.adjust(lrt$table$PValue, method="BH")
padj_lrt <- cbind(lrt$table,FDR)
write.csv(lrt$fitted.values,file = "edgeR_norm_value.csv",row.names = T)
write.csv(padj_lrt,file = "edgeR_results.csv",row.names = T)

FDR_0.05 = padj_lrt[padj_lrt$FDR < 0.05,]
FDR_0.05$up_down = ifelse(FDR_0.05$logFC>0,"up","down")
write.csv(FDR_0.05,file = "edgeR_FDR_0.05.csv",row.names = T)

up_gene <- FDR_0.05[FDR_0.05$logFC >= 1,]
down_gene <- FDR_0.05[FDR_0.05$logFC <= -1,]
dim(up_gene)
dim(down_gene)
write.csv(up_gene,file = "edgeR_up.csv",row.names = T)
write.csv(down_gene,file = "edgeR_down.csv",row.names = T)
