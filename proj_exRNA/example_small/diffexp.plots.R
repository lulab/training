## differential expression analysis
#-------------------------------
# filter expression matrix
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

