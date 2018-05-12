## sr2017_GSE71008 datasets

library(SingleCellExperiment)
library(scater)
options(stringsAsFactors = FALSE)

mx <- read.table("GSE71008.NvsCRC.reads.txt", sep = "\t")
anno <- read.table("GSE71008.NvsCRC.anno.txt", sep = "\t", header=T)

anno$Class <- "NC"
anno[which(anno$Stage=="1S"),]$Class <- "S1"
anno[which(anno$Stage=="2S"),]$Class <- "S2"
anno[which(anno$Stage=="3S"),]$Class <- "S3"
anno[which(anno$Stage=="4S"),]$Class <- "S4"

# construct singleCellExperiment object
anno_NvsEachStage <- anno
mx_NvsEachStage <- mx
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

# sample filtering

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

# save the data
saveRDS(reads_NvsEachStage, file = "GSE71008.reads_NvsEachStage.clean.rds")

# visualization
# log2 transfer
assay(reads_NvsEachStage, "logcounts_raw") <- log2(counts(reads_NvsEachStage) + 1)
reads_NvsEachStage.qc <- reads_NvsEachStage[rowData(reads_NvsEachStage)$use, colData(reads_NvsEachStage)$use]

# Before QC
endog_genes <- !rowData(reads_NvsEachStage)$is_feature_control
plotPCA(
    reads_NvsEachStage[endog_genes, ],
    exprs_values = "counts",
    colour_by = "Class",
    size_by = "total_features"
)

plotTSNE(
    reads_NvsEachStage[endog_genes, ],
    exprs_values = "counts",
    perplexity = 40,
    colour_by = "Class",
    size_by = "total_features",
    rand_seed = 123456,
    ntop = 100
)


# After QC
reads_NvsEachStage.qc <- reads_NvsEachStage[rowData(reads_NvsEachStage)$use, colData(reads_NvsEachStage)$use]
endog_genes <- !rowData(reads_NvsEachStage.qc)$is_feature_control
plotPCA(
    reads_NvsEachStage.qc[endog_genes, ],
    exprs_values = "counts",
    colour_by = "Class",
    size_by = "total_features"
)

plotTSNE(
    reads_NvsEachStage.qc[endog_genes, ],
    exprs_values = "counts",
    perplexity = 40,
    colour_by = "Class",
    size_by = "total_features",
    rand_seed = 123456,
    ntop = 100
)

# normalizaiton

# CPM
logcounts(reads_NvsEachStage.qc) <- log2(calculateCPM(reads_NvsEachStage.qc, use.size.factors = FALSE) + 1)
plotPCA(
    reads_NvsEachStage.qc[endog_genes, ],
    exprs_values = "logcounts",
    colour_by = "Class",
    size_by = "total_features"
)
plotRLE(
    reads_NvsEachStage.qc[endog_genes, ], 
    exprs_mats = list(Raw = "counts", CPM = "logcounts"),
    exprs_logged = c(TRUE, TRUE),
    colour_by = "Class"
)

# SF (DESeq)
reads_NvsEachStage.qc <- normaliseExprs(
    reads_NvsEachStage.qc,
    method = "RLE", 
    feature_set = endog_genes,
    return_log = TRUE,
    return_norm_as_exprs = TRUE
)
plotPCA(
    reads_NvsEachStage.qc[endog_genes, ],
    exprs_values = "normcounts",
    colour_by = "Class",
    size_by = "total_features"
)
plotRLE(
    reads_NvsEachStage.qc[endog_genes, ], 
    exprs_mats = list(Raw = "counts", SF = "normcounts"),
    exprs_logged = c(TRUE, TRUE),
    colour_by = "Class"
)

# TMM (edgeR)
reads_NvsEachStage.qc <- normaliseExprs(
    reads_NvsEachStage.qc,
    method = "TMM",
    feature_set = endog_genes,
    return_log = TRUE,
    return_norm_as_exprs = TRUE
)

plotPCA(
    reads_NvsEachStage.qc[endog_genes, ],
    exprs_values = "normcounts",
    colour_by = "Class",
    size_by = "total_features"
)

plotRLE(
    reads_NvsEachStage.qc[endog_genes, ], 
    exprs_mats = list(Raw = "counts", TMM = "normcounts"),
    exprs_logged = c(TRUE, TRUE),
    colour_by = "Class"
)


# scran
library(scran)
# define cluster for each sample
sampleLables <- c()
for(i in colnames(reads_NvsEachStage.qc)){tmp <- as.character(anno[which(anno$Individual==i),"Class"]); sampleLables <- c(sampleLables,tmp)}
sampleLables <- replace(sampleLables, which(sampleLables=="S1"),1)
sampleLables <- replace(sampleLables, which(sampleLables=="S2"),2)
sampleLables <- replace(sampleLables, which(sampleLables=="S3"),3)
sampleLables <- replace(sampleLables, which(sampleLables=="S4"),4)
sampleLables <- replace(sampleLables, which(sampleLables=="NC"),5)
sampleLables <- as.numeric(sampleLables)

# sampleLables <- quickCluster(reads_NvsEachStage.qc, min.size = 10)
reads_NvsEachStage.qc <- computeSumFactors(reads_NvsEachStage.qc, sizes = 10, clusters = sampleLables)
reads_NvsEachStage.qc <- normalize(reads_NvsEachStage.qc)

plotPCA(
    reads_NvsEachStage.qc, 
    exprs_values = "logcounts", 
    colour_by = "Class", 
    size_by = "total_features"
)
plotTSNE(
    reads_NvsEachStage.qc, 
    exprs_values = "logcounts", 
    perplexity = 10, 
    colour_by = "Class", 
    size_by = "total_features", 
    rand_seed = 123456)
plotRLE(
    reads_NvsEachStage.qc, 
    exprs_mats = list(Raw = "counts", scran = "logcounts"), 
    exprs_logged = c(TRUE, TRUE), 
    colour_by = "Class"
)


# identify confounding factors
# TBD

