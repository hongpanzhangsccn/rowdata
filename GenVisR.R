#source("https://bioconductor.org/biocLite.R")
#biocLite("GenomeInfoDbData")
#biocLite("GenVisR")

setwd("C:\\Users\\lexb4\\Desktop\\SNP\\05.waterfall")

library(GenVisR)
rt=read.table("waterfallInput.txt",header=T,sep="\t",check.names=F,quote="")
pdf(file="waterfall.pdf",height=9,width=12)
waterfall(rt)
dev.off()



