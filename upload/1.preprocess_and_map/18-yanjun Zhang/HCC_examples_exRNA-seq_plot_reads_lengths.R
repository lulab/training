Args <- commandArgs(TRUE)
infile <- Args[1]
outfile <- Args[2]

mx <- read.table(infile,header = T,sep="\t")

library(ggplot2)
p <- ggplot(mx, aes(len,num)) + geom_col() + theme_bw() + facet_grid(type ~.)
pdf(outfile,width=4,height=8)
plot(p)
dev.off()
