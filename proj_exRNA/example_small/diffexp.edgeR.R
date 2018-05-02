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

#--------------------
# edgeR
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
write.table(tidyResult,file="hcc_example.miRNA.NCvsHCC.edgeR.classic.tsv",sep="\t",row.names=FALSE)

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


