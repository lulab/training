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
# DESeq2:
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

# generate Dataset 
dds <- DESeqDataSetFromMatrix(countData, colData, design=~Treatment, tidy=TRUE)

# normlize using rlog mathod
norm <- rlog(dds,blind=FALSE)
norm_matrix <- assay(norm)
norm_df <- data.frame(Gene=rownames(norm_matrix), norm_matrix)
write.table(norm_df, "hcc_example.miRNA.homer.DESeq2.rlog.mx", row.names = FALSE,sep="\t")

deg <- DESeq(dds)
res <- results(deg,tidy=TRUE)
merged_res <- merge(norm_df,res,by.x="Gene",by.y="row")
write.table(merged_res,file="hcc_example.miRNA.NCvsHCC.DESeq2.tsv",sep="\t",row.names=FALSE)










