# R script to find differentially expressed genes by R package wilcoxon.

library("limma")
library("edgeR")

d<-read.csv("CPM_filtered_outlier.csv",header = T)
RPM <- d[,2:42]
rownames(RPM)<-d[,1]
keep <- rowSums(RPM>0)>=9
RPM <- RPM[keep,]
group <- c(rep("HCC",30),rep("Normal",11))

RPM_HCC <- RPM[,1:30]
RPM_Normal <- RPM[,31:41]
RPM_HCC_avg <- apply(RPM_HCC,1,mean)
RPM_Normal_avg <- apply(RPM_Normal,1,mean)
RPM_Avg <- cbind(RPM_HCC_avg,RPM_Normal_avg)
write.csv(RPM_Avg,file = "RPM_average.csv",row.names = T)

myFun <- function(x){
  x = as.vector(x)
  v1 = x[1:30]
  v2 = x[31:41]
  out <- wilcox.test(v1,v2)
  out <- out$p.value
}
p_value <- apply(RPM,1,myFun)

d2<-read.csv("RPM_average_log2FC.csv",header = T)
Result <- d2[,2:4]
Result<- cbind(Result,p_value)
FDR <- p.adjust(Result[,4],method = "BH")
Result <- cbind(Result,FDR)
write.csv(Result,file = "wilcoxon_result.csv",row.names = T)

FDR_0.05 = Result[Result$FDR<0.05,]
FDR_0.05$up_down = ifelse(as.numeric(FDR_0.05$log2FoldChange) > 0,"up","down")
write.csv(FDR_0.05,file = "wilcoxon_FDR_0.05.csv",row.names = T)

up_gene <- FDR_0.05[as.numeric(FDR_0.05$log2FoldChange) >= 1,]
down_gene <- FDR_0.05[as.numeric(FDR_0.05$log2FoldChange) <= -1,]
dim(up_gene)
dim(down_gene)
write.csv(up_gene,file = "wilcoxon_up.csv",row.names = T)
write.csv(down_gene,file = "wilcoxon_down.csv",row.names = T)
