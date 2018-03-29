# R script to calculate gene's cpm and create a heatmap nny heatmap.2
```
counts <- read.table("HCC_ex_miRNA_primary.counts",header = T, sep = "\t")
d<-counts[,4:23]
rownames(d)<-counts[,2]
group<-factor(c("B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","A","A","A","A","N"))
design <-model.matrix(~group)
d<-d[keep,,keep.lib.sizes=F]
CPM<-cpm(d)
write.csv(CPM,file="HCC_ex_miRNA_primary.cpm")
CPM_matrix <-as.matrix(CPM)
sample_color <-c(rep("#FF5733",15),rep("#FFC300",4),"#74C753")
gene_color <-c(rep("#86c0e8",46),rep("#DAE4EB",14))
heatmap.2(CPM, col=colorpanel(100,low = "white",high = "steelblue"), key = T, symkey = F, density.info = "none", trace = "none", cexRow = 1, cexCol = 1, srtCol = 45, ColSideColors = sample_color, RowSideColors = gene_color, scale = "row", keysize = 1, margins = c(5,10), dendrogram = c("both"))
```
