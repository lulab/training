#!/bin/bash

path0=/home/younglee/projects/hcc_example
mkdir $path0/05.diffexp/

## differential expression analysis
#-------------------------------
# 0.filter expression matrix
mx <- read.table("hcc_example.miRNA.homer.ct.mx",sep="\t",header=T)
head(mx[,1:5])

# check library size (total mapped reads)
colSums(mx[,2:ncol(mx)])

# check the number of detected genes
apply(mx[,2:ncol(mx)], 2, function(c)sum(c!=0))

# filter genes
filter_genes <- apply(
    mx[,2:ncol(mx)],
    1,
    function(x) length(x[x > 2]) >= 2
)

mx_filterGenes <- mx[filter_genes,]
head(mx_filterGenes[,1:5])

#--------------------------
# 1.DESeq2:
# experimential design
design <- read.table("design.tsv",sep="\t",header=T)

# Sample  Treatment
# NC_1    NC
# NC_2    NC
# NC_3    NC
# BeforeSurgery_1 HCC
# BeforeSurgery_2 HCC
# BeforeSurgery_3 HCC
# AfterSurgery_1  HCC
# AfterSurgery_2  HCC
# AfterSurgery_3  HCC


# expression matrix
head(mx_filterGenes[,1:5])
#     geneID NC_1 NC_2 NC_3 BeforeSurgery_1
#  MI0016824    0    0    0               3
#  MI0017385    9    5   10              16
#  MI0016429    6    8    1               1
#  MI0003676 1021  197  713              35
#  MI0025922    3    1    3               8
#  MI0016849    1    1    2               3


# basic script for normalizing with DESeq2
library(DESeq2)
#Read Data in
countData <- mx_filterGenes
colData <- design

# generate Dataset object 
dds <- DESeqDataSetFromMatrix(countData, colData, design=~Treatment, tidy=TRUE)

# normlize using rlog mathod
norm <- rlog(dds,blind=FALSE)
norm_matrix <- assay(norm)
norm_df <- data.frame(Gene=rownames(norm_matrix), norm_matrix)
write.table(norm_df, "hcc_example.miRNA.homer.DESeq2.rlog.mx", quote=F, row.names = FALSE,sep="\t")

deg <- DESeq(dds)
res <- results(deg,tidy=TRUE)
merged_res <- merge(norm_df,res,by.x="Gene",by.y="row")
write.table(merged_res,file="hcc_example.miRNA.NCvsHCC.DESeq2.tsv",quote=F, sep="\t",row.names=FALSE)

# MA-plot
pdf("hcc_example.miRNA.NCvsHCC.DESeq2.MAplot.pdf")
plotMA(deg, ylim=c(-5,5))
dev.off()

# visualize one diff-exp gene
pdf("hcc_example.miRNA.NCvsHCC.DESeq2.plotCounts.pdf")
plotCounts(deg, gene=which.min(res$padj), intgroup="Treatment")
dev.off()

# heatmap for random genes
pdf("hcc_example.miRNA.NCvsHCC.DESeq2.pheatmap.pdf")
library("pheatmap")
ntd <- normTransform(deg)
select <- order(rowMeans(counts(deg,normalized=TRUE)), decreasing=TRUE)[1:25]
df <- as.data.frame(colData(deg)[,c("Treatment","Sample")])
pheatmap(assay(norm)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=TRUE, annotation_col=df)
dev.off()

# pca 
pdf("hcc_example.miRNA.NCvsHCC.DESeq2.pca.pdf")
plotPCA(norm, intgroup=c("Treatment"))
dev.off()

# distance between samples
pdf("hcc_example.miRNA.NCvsHCC.DESeq2.sampleDists.pdf")
sampleDists <- dist(t(assay(norm)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(norm$Treatment, norm$Sample, sep="-")
colnames(sampleDistMatrix) <- paste(norm$Treatment, norm$Sample, sep="-")
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()


#--------------------
# 2.edgeR
# basic script for running edgeR
library(edgeR)
#Read Data in
countData <- mx_filterGenes[,-1]
rownames(countData) <- mx_filterGenes[,1]
design <- read.table("design.tsv",sep="\t",header=T)
colData <- design

# generate DGE object
y <- DGEList(countData, samples=colData, group=colData$Treatment)
y <- calcNormFactors(y)

#Estimate Error Model
design <-model.matrix(~Treatment, data=colData)
y <- estimateDisp(y, design)

# classic methods: compute p-values, then output
et <- exactTest(y)
res <- topTags(et,Inf)
tidyResult <- data.frame(Gene=rownames(res$table), res$table)
write.table(tidyResult,file="hcc_example.miRNA.NCvsHCC.edgeR.classic.tsv", quote=F, sep="\t",row.names=FALSE)

# Generalized linear models 
fit <- glmFit(y,design)
# likelihood ratio test
lrt <- glmLRT(fit,contrast = c(1,-1))
FDR <- p.adjust(lrt$table$PValue, method="BH")
padj_lrt <- cbind(lrt$table,FDR)
fit_df <- lrt$fitted.values
write.table(fit_df,file = "hcc_example.miRNA.homer.edgeR.TMM.mx",row.names = T, sep="\t", quote=F)
merged_lrt <- merge(fit_df,padj_lrt,by="row.names")
colnames(merged_lrt)[1] <- "Genes"
write.table(merged_lrt,file = "hcc_example.miRNA.NCvsHCC.edgeR.tsv",row.names = F, sep="\t", quote=F)

#----------------
# 3.using homer
getDiffExpression.pl ../04.counts/hcc_example.miRNA.homer.ct.tsv NC NC NC HCC HCC HCC HCC HCC HCC -repeats -DESeq2 > hcc_example.miRNA.NCvsHCC.homer.DESeq2.tsv
getDiffExpression.pl ../04.counts/hcc_example.miRNA.homer.ct.tsv NC NC NC HCC HCC HCC HCC HCC HCC -repeats -edgeR > hcc_example.miRNA.NCvsHCC.homer.edgeR.tsv


# 4. limma
library(limma)
countData <- mx_filterGenes[,-1]
rownames(countData) <- mx_filterGenes[,1]
design <- model.matrix(~ 0+factor(c(1,1,1,2,2,2,2,2,2)))
colnames(design) <- c("NC", "HCC")

# generate DGE object
dge <- DGEList(countData)
keep <- filterByExpr(dge, design)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)

# limma-trend model
# If the sequencing depth is reasonably consistent across the RNA samples, then the simplest and most robust approach to differential exis to use limma-trend.
logCPM <- cpm(dge, log=TRUE, prior.count=0.25)
fit <- lmFit(logCPM, design)
fit <- eBayes(fit, trend=TRUE)
topTable(fit, coef=ncol(design))

# voom model
# the library sizes are quite variable between samples, then the voom approach is theoretically more powerful than limma-trend.
v <- voom(dge, design, plot=TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit)
topTable(fit, coef=ncol(design))

# 5. wilcox test 
cpmMx <- read.table("hcc_example.miRNA.homer.rpm.mx",sep="\t",header=T)
filter_cpm <- apply(
    mx[,2:ncol(cpmMx)],
    1,
    function(x) length(x[x > 0]) >= 2
)
mx_filterCPM <- cpmMx[filter_cpm,]

myFun <- function(x){
  x = as.numeric(x)
  v1 = x[2:4]
  v2 = x[5:10]
  out <- wilcox.test(v1,v2)
  out <- out$p.value
}
p_value <- apply(mx_filterCPM,1,myFun)
p_value[is.nan(p_value)] <- 1
FDR <- p.adjust(p_value,method = "BH")
mx_filterCPM$avgNC <- apply(foo[,2:4],1,mean)
mx_filterCPM$avgHCC <- apply(foo[,5:10],1,mean)
mx_filterCPM$log2fc <- log2((foo$avgNC+0.25)/(foo$avgHCC+0.25))
results <- cbind(mx_filterCPM,p_value,FDR)
write.table(results,file = "hcc_example.miRNA.NCvsHCC.wilcox.tsv",row.names = F, sep="\t", quote=F)

#--------------------------------------------
## filter out differential expressed genes

awk 'BEGIN{FS=OFS="\t"}($12>1 && $15<=0.05){print $1}' hcc_example.miRNA.NCvsHCC.DESeq2.tsv > hcc_example.miRNA.NCvsHCC.DESeq2.NC_high.ids
awk 'BEGIN{FS=OFS="\t"}($12<-1 && $15<=0.05){print $1}' hcc_example.miRNA.NCvsHCC.DESeq2.tsv > hcc_example.miRNA.NCvsHCC.DESeq2.HCC_high.ids

mx <- read.table("hcc_example.miRNA.NCvsHCC.DESeq2.ids.rpkm.forPlot",sep="\t",head=T)
# log2 transfer the RPKM value
log2mx <- log2(mx)
library(gplots)
pdf("hcc_example.miRNA.NCvsHCC.DESeq2.ids.rpkm.heatmap.pdf",width=4,height=6)
heatmap.2(as.matrix(log2mx), Rowv=T, Colv=T, dendrogram = "both", trace = "none", density.info = c("none"), scale="row")
dev.off()




