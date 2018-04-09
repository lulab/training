#This Rscript was written for normalizing raw counts to cpm and choosing differentially expressed genes between two groups accodring to pvalues of wilcoxon test 

Args=commandArgs(TRUE)

raw_counts <-read.table(Args[1], header=T, sep='\t')
depth <- read.table(Args[2] ,head=F, sep='\t')
norm_counts <- raw_counts
for (i in 1:dim(depth)[1]) {norm_counts[,i]=raw_counts[,i]/depth[i,1]*1000000}

c<-data.frame(t(norm_counts),group = factor(c(rep(Args[3],each=Args[4]),rep(Args[5],each=Args[6]))))

out=c();
for(i in 1:(dim(c)[2]-1)){
    tmp=wilcox.test(c[,i] ~ group,paired = FALSE, c)$p.value;
    out=c(out,tmp);
}
norm_counts$wilcox <- out
hi <- norm_counts[norm_counts$wilcox <0.05,]
write.table(hi,Args[7], quote=FALSE,sep='\t')
