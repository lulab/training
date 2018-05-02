## differential expression analysis
#-------------------------------
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





