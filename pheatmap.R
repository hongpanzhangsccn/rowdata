#install.packages("pheatmap")

setwd("C:\\Users\\lexb4\\Desktop\\TCGAimmune\\09.heatmap")      #设置工作目录
rt=read.table("CIBERSORT.filter.txt",sep="\t",header=T,row.names=1,check.names=F)
rt=t(rt)

library(pheatmap)
Type=c(rep("con",4),rep("treat",124))    #修改对照和处理组样品数目
names(Type)=colnames(rt)
Type=as.data.frame(Type)

pdf("heatmap.pdf",height=6,width=15)
pheatmap(rt, 
         annotation=Type, 
         color = colorRampPalette(c("green", "black", "red"))(50),
         cluster_cols =F,
         fontsize = 8,
         fontsize_row=7,
         fontsize_col=6)
dev.off()

