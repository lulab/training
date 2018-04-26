# R script
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

# check the correlations between each samples
cor.test(mx_filterGenes$NC_1,mx_filterGenes$NC_2,methods="spearman")

pdf("NC_cor.pdf")
pairs(~NC_1+NC_2+NC_3,data=mx_filterGenes)
dev.off()

# save the results
write.table(mx_filterGenes,"hcc_example.miRNA.homer.ct.filtered.mx",sep="\t",quote=F,col.names=T,row.names=F)
