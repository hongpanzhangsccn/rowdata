#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("karyoploteR")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("pasillaBamSubset")

#设置工作目录
setwd("C:\\Users\\lexb4\\Desktop\\ICGCsnp\\08.mutationPos")
library(karyoploteR)
rt=read.table("mutStatPos.txt",sep="\t",header=T)
rt[,3]=rt[,3]+1000000
rt=rt[rt$Num>5,]
dd=makeGRangesFromDataFrame(rt)
pdf(file="mutationPos.pdf",width=10,height=7)
kp <- plotKaryotype("hg19", plot.type=1)
y1=log2(rt[,4]+1)
kpHeatmap(kp, dd, y=y1, colors = c("green", "black", "red"), r0=0, r1=0.65)
kpAddBaseNumbers(kp, tick.dist=10000000, minor.tick.dist=1000000)
dev.off()
