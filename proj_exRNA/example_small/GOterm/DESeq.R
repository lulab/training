# filter expression matrix
mx <- read.table("hcc_example.mRNA_exon.homer.ct.mx",sep="\t",header=T)

# filter genes
filter_genes <- apply(
    mx[,2:ncol(mx)],
    1,
    function(x) length(x[x > 2]) >= 2
)

mx_filterGenes <- mx[filter_genes,]

# load experimential design
design <- read.table("design.tsv",sep="\t",header=T)

#-----------------------------------------
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
write.table(norm_df, "hcc_example.mRNA_exon.homer.DESeq2.rlog.mx", row.names = FALSE,sep="\t")

deg <- DESeq(dds)
res <- results(deg,tidy=TRUE)
merged_res <- merge(norm_df,res,by.x="Gene",by.y="row")
write.table(merged_res,file="hcc_example.mRNA_exon.NCvsHCC.DESeq2.tsv",sep="\t",row.names=FALSE)
