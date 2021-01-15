#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("edgeR")

#install.packages("pheatmap")


#引用包
library(edgeR)
library(pheatmap)

logFCfiler=2                  #logFC过滤条件
fdrFilter=0.05                #矫正后的p值过滤条件
conNum=353                    #野生型样品数目
treatNum=22                   #突变型样品数目

setwd("C:\\Users\\lexb4\\Desktop\\snpRNA\\11.diff")                    #设置工作目录
rt=read.table("sampleExp.txt",sep="\t",header=T,check.names=F)         #读取输入文件

#如果一个基因出现多行，取均值
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>1,]

#对数据进行矫正，差异分析
group=c(rep("normal",conNum),rep("tumor",treatNum))
design <- model.matrix(~group)
y <- DGEList(counts=data,group=group)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y,pair = c("normal","tumor"))
ordered_tags <- topTags(et, n=100000)

allDiff=ordered_tags$table
allDiff=allDiff[is.na(allDiff$FDR)==FALSE,]
diff=allDiff
newData=y$pseudo.counts

#输出差异结果
write.table(diff,file="all.xls",sep="\t",quote=F)
diffSig = diff[(diff$FDR < fdrFilter & (diff$logFC>logFCfiler | diff$logFC<(-logFCfiler))),]
diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut, file="diff.xls",sep="\t",quote=F,col.names=F)
write.table(diffSigOut, file="diff.txt",sep="\t",quote=F,col.names=F)

#输出矫正后的文件
normalizeExp=rbind(id=colnames(newData),newData)
write.table(normalizeExp,file="normalExp.txt",sep="\t",quote=F,col.names=F)   #输出所有基因校正后的表达值
diffExp=rbind(id=colnames(newData),newData[rownames(diffSig),])
write.table(diffExp,file="diffExp.txt",sep="\t",quote=F,col.names=F)         #输出差异基因校正后的表达值

#绘制差异基因热图
geneNum=20
diffSig=diffSig[order(as.numeric(as.vector(diffSig$logFC))),]
diffGeneName=as.vector(rownames(diffSig))
diffLength=length(diffGeneName)
hmGene=c()
if(diffLength>40){
    hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
}else{
    hmGene=diffGeneName
}
hmExp=newData[hmGene,]
hmExp=log2(hmExp+0.001)
Type=c(rep("Wild",conNum),rep("Mut",treatNum))
names(Type)=colnames(newData)
Type=as.data.frame(Type)
pdf(file="heatmap.pdf",height=4,width=7)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("green", "black", "red"))(50),
         cluster_cols =F,
         show_colnames = F,
         #scale="row",
         fontsize = 8,
         fontsize_row=6,
         fontsize_col=8)
dev.off()

#绘制火山图
pdf(file="vol.pdf",width = 5,height =5)
xMax=10
yMax=20
plot(allDiff$logFC,-log10(allDiff$FDR), ylab="-log10(FDR)",xlab="logFC",
     main="Volcano", ylim=c(0,yMax),xlim=c(-xMax,xMax),yaxs="i",pch=20, cex=1)
diffSub=allDiff[allDiff$FDR<fdrFilter & allDiff$logFC>logFCfiler,]
points(diffSub$logFC, -log10(diffSub$FDR), pch=20, col="red",cex=1)
diffSub=allDiff[allDiff$FDR<fdrFilter & allDiff$logFC<(-logFCfiler),]
points(diffSub$logFC, -log10(diffSub$FDR), pch=20, col="green",cex=1)
abline(v=0,lty=2,lwd=3)
dev.off()


