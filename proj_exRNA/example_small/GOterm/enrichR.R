# using enrichR to preform GO term analysis
Args <- commandArgs(TRUE);
infile <- Args[1];
outprefix <- Args[2];

mx <- read.table(infile,sep="\t",header=F)
genelist <- as.vector(mx$V1)
rm(mx)

library(enrichR)
dbs <- listEnrichrDbs()
dbs <- c("GO_Molecular_Function_2017b", "GO_Cellular_Component_2017b", "GO_Biological_Process_2017b", "KEGG_2016")
enriched <- enrichr(genelist, dbs)
printEnrich(enriched, paste(outprefix,".enrichR.txt",sep=""), sep = "\t", columns = c(1:9))

bp_out <- enriched[["GO_Biological_Process_2017b"]]
bp_out$log10P <- -log10(bp_out$Adjusted.P.value)
cc_out <- enriched[["GO_Cellular_Component_2017b"]]
cc_out$log10P <- -log10(cc_out$Adjusted.P.value)
mf_out <- enriched[["GO_Molecular_Function_2017b"]]
mf_out$log10P <- -log10(mf_out$Adjusted.P.value)
kegg_out <- enriched[["KEGG_2016"]]
kegg_out$log10P <- -log10(kegg_out$Adjusted.P.value)

write.table(bp_out,paste(outprefix,".enrichR.bp.txt", sep=""),sep="\t",quote=F,col.names=T,row.names=F)
write.table(cc_out,paste(outprefix,".enrichR.cc.txt", sep=""),sep="\t",quote=F,col.names=T,row.names=F)
write.table(mf_out,paste(outprefix,".enrichR.mf.txt", sep=""),sep="\t",quote=F,col.names=T,row.names=F)
write.table(kegg_out,paste(outprefix,".enrichR.kegg.txt", sep=""),sep="\t",quote=F,col.names=T,row.names=F)

bp_selt <- bp_out[order(bp_out$log10P, decreasing=T)[1:20],]
bp_selt <- subset(bp_selt,bp_selt$Adjusted.P.value<=0.05)
if(nrow(bp_selt)>=1){bp_selt$class <- "Biological_Process"}
cc_selt <- cc_out[order(cc_out$log10P, decreasing=T)[1:20],]
cc_selt <- subset(cc_selt,cc_selt$Adjusted.P.value<=0.05)
if(nrow(cc_selt)>=1){cc_selt$class <- "Cellular_Component"}
mf_selt <- mf_out[order(mf_out$log10P, decreasing=T)[1:20],]
mf_selt <- subset(mf_selt,mf_selt$Adjusted.P.value<=0.05)
if(nrow(mf_selt)>=1){mf_selt$class <- "Molecular_Function"}
kegg_selt <- kegg_out[order(kegg_out$log10P, decreasing=T)[1:20],]
kegg_selt <- subset(kegg_selt,kegg_selt$Adjusted.P.value<=0.05)
if(nrow(kegg_selt)>=1){kegg_selt$class <- "KEGG_Pathway"}

mainTitle <- tail(unlist(strsplit(outprefix, "/")),n=1)
merge_out <- rbind(bp_selt,cc_selt,mf_selt,kegg_selt)
merge_out$regions <- mainTitle
write.table(merge_out,paste(outprefix,".enrichR.selt.txt", sep=""),sep="\t",quote=F,col.names=T,row.names=F)

library(ggplot2)
# Basic barplot
pdf(paste(outprefix,".enrichR.barplot.pdf",sep=""),width=8,height=6)
p <- ggplot(data=merge_out, aes(x=Term, y=log10P)) + geom_bar(stat="identity") + coord_flip() + theme_bw() + theme(axis.text.x = element_text(size=5)) + facet_grid(class ~ ., scales = "free", space = "free") + labs(title = mainTitle)
plot(p)
dev.off()


