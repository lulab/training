## sr2017_GSE71008 datasets

library(SingleCellExperiment)
library(scater)
options(stringsAsFactors = FALSE)
setwd("./GSE71008")

mx <- read.table("GSE71008.NvsCRC.reads.txt", sep = "\t")
anno <- read.table("GSE71008.NvsCRC.annoWithBatch.txt", sep = "\t", header=T)

# construct singleCellExperiment object
anno_NvsEachStage <- anno
mx_NvsEachStage <- mx[,colnames(mx) %in% anno$Individual]
reads_NvsEachStage <- SingleCellExperiment(
    assays = list(counts = as.matrix(mx_NvsEachStage)),
    colData = anno_NvsEachStage)

# filter samples
keep_feature <- rowSums(counts(reads_NvsEachStage) > 0) > 0
reads_NvsEachStage <- reads_NvsEachStage[keep_feature, ]

stableRNA <- isSpike(reads_NvsEachStage, "stableRNA") <- rownames(reads_NvsEachStage) %in% c("mature_miRNA:hsa-miR-99a-5p", "mature_miRNA:hsa-miR-30a-5p", "mature_miRNA:hsa-miR-221-3p")

reads_NvsEachStage <- calculateQCMetrics(
    reads_NvsEachStage,
    feature_controls = list(
        stableRNA = isSpike(reads_NvsEachStage, "stableRNA")
    )
)

hist(reads_NvsEachStage$total_counts,breaks = 100)
abline(v=990000, col="red")
filter_by_total_counts <- (reads_NvsEachStage$total_counts > 990000)
table(filter_by_total_counts)

hist(reads_NvsEachStage$total_features,breaks = 100)
abline(v=2500, col="red")
filter_by_expr_features <- (reads_NvsEachStage$total_features > 2500)
table(filter_by_expr_features)


plotPhenoData(
    reads_NvsEachStage,
    aes_string(
        x = "total_features",
        y = "pct_counts_stableRNA",
        colour = "Class"
    )
)

filter_by_endoCtrl <- reads_NvsEachStage$pct_counts_stableRNA < 10
table(filter_by_endoCtrl)

reads_NvsEachStage$use <- (
    # sufficient features (genes)
    filter_by_expr_features &
    # sufficient molecules counted
    filter_by_total_counts &
    # sufficient endogenous RNA
    filter_by_endoCtrl
)
table(reads_NvsEachStage$use)

## gene QC
# filter genes with too low expression
plotQC(reads_NvsEachStage, type = "highest-expression")

filter_genes <- apply(counts(reads_NvsEachStage[, colData(reads_NvsEachStage)$use]), 1, function(x) length(x[x >= 2]) >= 20)
table(filter_genes)
rowData(reads_NvsEachStage)$use <- filter_genes
dim(reads_NvsEachStage[rowData(reads_NvsEachStage)$use, colData(reads_NvsEachStage)$use])

# visualization
# log2 transfer
assay(reads_NvsEachStage, "logcounts_raw") <- log2(counts(reads_NvsEachStage) + 1)
reads_NvsEachStage.qc <- reads_NvsEachStage[rowData(reads_NvsEachStage)$use, colData(reads_NvsEachStage)$use]

# save the data
saveRDS(reads_NvsEachStage.qc, file = "GSE71008.reads_NvsEachStage.clean.batch.rds")


##---------------------
## Imputation
library("scImpute")
library(scran)
reads_NvsEachStage.qc <- readRDS("GSE71008.reads_NvsEachStage.clean.batch.rds")
anno <- read.table("GSE71008.NvsCRC.annoWithBatch.txt", sep = "\t", header=T)

write.csv(counts(reads_NvsEachStage.qc), "GSE71008.reads_NvsEachStage.qc.rds.csv")

scimpute(count_path = "GSE71008.reads_NvsEachStage.qc.rds.csv", infile = "csv", outfile = "txt", out_dir = "./", ncores = 2, Kcluster = 5)
res.qc <- read.table("./scimpute_count.txt")
reads_NvsEachStage.qc.impute <- SingleCellExperiment(assays = list(counts = as.matrix(res.qc)), colData = colData(reads_NvsEachStage.qc))
reads_NvsEachStage.qc.impute <- calculateQCMetrics(reads_NvsEachStage.qc.impute)
assay(reads_NvsEachStage.qc.impute, "logcounts_raw") <- log2(counts(reads_NvsEachStage.qc.impute) + 1)

# compare the raw and after imputation

plotPCA(
    reads_NvsEachStage.qc,
    exprs_values = "counts",
    colour_by = "Class",
    size_by = "total_features"
) + ggtitle("reads_NvsEachStage.qc")

plotPCA(
    reads_NvsEachStage.qc.impute,
    exprs_values = "counts",
    colour_by = "Class",
    size_by = "total_features"
) + ggtitle("reads_NvsEachStage.qc.impute")


plotTSNE(
    reads_NvsEachStage.qc,
    exprs_values = "counts",
    perplexity = 10,
    colour_by = "Class",
    size_by = "total_features",
    rand_seed = 123456
)

plotTSNE(
    reads_NvsEachStage.qc.impute,
    exprs_values = "counts",
    perplexity = 10,
    colour_by = "Class",
    size_by = "total_features",
    rand_seed = 123456
)


# scran
library(scran)

# Cluster similar cells based on rank correlations in their gene expression profiles
qclust <- quickCluster(reads_NvsEachStage.qc.impute, min.size = 10)
# Methods to normalize single-cell RNA-seq data by deconvolving size factors from cell pools.
reads_NvsEachStage.qc.impute <- computeSumFactors(reads_NvsEachStage.qc.impute, sizes = 10, clusters = qclust)
reads_NvsEachStage.qc.impute <- normalize(reads_NvsEachStage.qc.impute)

plotPCA(
    reads_NvsEachStage.qc.impute,
    exprs_values = "logcounts",
    colour_by = "Class",
    size_by = "total_features"
) + ggtitle("reads_NvsEachStage.qc.impute.scran")

plotTSNE(
    reads_NvsEachStage.qc.impute,
    exprs_values = "logcounts",
    perplexity = 10,
    colour_by = "Class",
    size_by = "total_features",
    rand_seed = 123456
)

plotRLE(
    reads_NvsEachStage.qc.impute,
    exprs_mats = list(Raw = "counts", SF = "logcounts"),
    exprs_logged = c(TRUE, TRUE),
    colour_by = "Class"
)

# SF (DESeq)
reads_NvsEachStage.qc.impute <- normaliseExprs(
    reads_NvsEachStage.qc.impute,
    method = "RLE", 
    return_log = TRUE,
    return_norm_as_exprs = TRUE
)
plotPCA(
    reads_NvsEachStage.qc.impute,
    exprs_values = "normcounts",
    colour_by = "Class",
    size_by = "total_features"
)
plotRLE(
    reads_NvsEachStage.qc.impute, 
    exprs_mats = list(Raw = "counts", SF = "normcounts"),
    exprs_logged = c(TRUE, TRUE),
    colour_by = "Class"
)

# TMM (edgeR)
reads_NvsEachStage.qc.impute <- normaliseExprs(
    reads_NvsEachStage.qc.impute,
    method = "TMM",
    return_log = TRUE,
    return_norm_as_exprs = TRUE
)

plotPCA(
    reads_NvsEachStage.qc.impute,
    exprs_values = "normcounts",
    colour_by = "Class",
    size_by = "total_features"
)

plotRLE(
    reads_NvsEachStage.qc.impute, 
    exprs_mats = list(Raw = "counts", TMM = "normcounts"),
    exprs_logged = c(TRUE, TRUE),
    colour_by = "Class"
)

# identify confounding factors
## correlations with PCs
plotPCA(
    reads_NvsEachStage.qc.impute,
    exprs_values = "counts",
    colour_by = "RNA_Isolation_batch",
    size_by = "total_features"
) + ggtitle("reads_NvsEachStage.qc.impute")

plotQC(
    reads_NvsEachStage.qc.impute,
    type = "find-pcs",
    exprs_values = "logcounts",
    variable = "RNA_Isolation_batch"
)
 
plotQC(
    reads_NvsEachStage.qc.impute,
    type = "expl",
    exprs_values = "logcounts",
    variables = c(
        "total_features",
        "total_counts",
        "RNA_Isolation_batch",
        "Class",
        "library_preparation_day"
    )
)

#mes(reads_NvsEachStage.qc.impute)# remove batch effect using enfogenous genes
library(RUVSeq)
library(sva)
library(scRNA.seq.funcs)

# using RUVseq (RUVs)
# RUVs uses centered (technical) replicate/negative control samples for which the covariates of interest are constant
scIdx <- matrix(-1, ncol = max(table(reads_NvsEachStage.qc.impute$Class)), nrow = 5)
tmp <- which(reads_NvsEachStage.qc.impute$Class == "S1")
scIdx[1, 1:length(tmp)] <- tmp
tmp <- which(reads_NvsEachStage.qc.impute$Class == "S2")
scIdx[2, 1:length(tmp)] <- tmp
tmp <- which(reads_NvsEachStage.qc.impute$Class == "S3")
scIdx[3, 1:length(tmp)] <- tmp
tmp <- which(reads_NvsEachStage.qc.impute$Class == "S4")
scIdx[4, 1:length(tmp)] <- tmp
tmp <- which(reads_NvsEachStage.qc.impute$Class == "NC")
scIdx[5, 1:length(tmp)] <- tmp


cIdx <- rownames(reads_NvsEachStage.qc.impute)
ruvs <- RUVs(logcounts(reads_NvsEachStage.qc.impute), cIdx, k = 1, scIdx = scIdx, isLog = TRUE)
assay(reads_NvsEachStage.qc.impute, "ruvs1") <- ruvs$normalizedCounts
ruvs <- RUVs(logcounts(reads_NvsEachStage.qc.impute), cIdx, k = 5, scIdx = scIdx, isLog = TRUE)
assay(reads_NvsEachStage.qc.impute, "ruvs5") <- ruvs$normalizedCounts
ruvs <- RUVs(logcounts(reads_NvsEachStage.qc.impute), cIdx, k = 10, scIdx = scIdx, isLog = TRUE)
assay(reads_NvsEachStage.qc.impute, "ruvs10") <- ruvs$normalizedCounts
ruvs <- RUVs(logcounts(reads_NvsEachStage.qc.impute), cIdx, k = 20, scIdx = scIdx, isLog = TRUE)
assay(reads_NvsEachStage.qc.impute, "ruvs20") <- ruvs$normalizedCounts
ruvs <- RUVs(logcounts(reads_NvsEachStage.qc.impute), cIdx, k = 25, scIdx = scIdx, isLog = TRUE)
assay(reads_NvsEachStage.qc.impute, "ruvs25") <- ruvs$normalizedCounts



pdf("GSE71008.reads_NvsEachStage.remove_confounders.pca.pdf")
for(n in assayNames(reads_NvsEachStage.qc.impute)) {
    print(
        plotPCA(
            reads_NvsEachStage.qc.impute,
            colour_by = "Class",
            size_by = "total_features",
            exprs_values = n
        ) +
        ggtitle(n)
    )
}
dev.off()


# Effectiveness 1: PCA
pdf("GSE71008.reads_NvsEachStage.remove_confounders.pca_Batch.pdf")
for(n in assayNames(reads_NvsEachStage.qc.impute)) {
    print(
        plotPCA(
            reads_NvsEachStage.qc.impute,
            colour_by = "RNA_Isolation_batch",
            size_by = "total_features",
            shape_by = "library_preparation_day",
            exprs_values = n
        ) +
        ggtitle(n)
    )
}
dev.off()

# Effectiveness 2: RLE
pdf("GSE71008.reads_NvsEachStage.remove_confounders.plotRLE.pdf")
res <- list()
for(n in assayNames(reads_NvsEachStage.qc.impute)) {
    res[[n]] <- suppressWarnings(calc_cell_RLE(assay(reads_NvsEachStage.qc.impute, n)))
}
par(mar=c(6,4,1,1))
boxplot(res, las=2)
dev.off()

# Effectiveness 3: plotQC
pdf("GSE71008.reads_NvsEachStage.remove_confounders.plotQC.pdf")
for(n in assayNames(reads_NvsEachStage.qc.impute)) {
    print(
        plotQC(
            reads_NvsEachStage.qc.impute,
            type = "expl",
            exprs_values = n,
            variables = c(
                "total_features",
                "total_counts",
                "RNA_Isolation_batch",
                "Class",
                "library_preparation_day"
            )
        ) +
        ggtitle(n)
    )
}
dev.off()



