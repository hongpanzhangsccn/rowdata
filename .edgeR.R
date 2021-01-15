#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("edgeR")

#install.packages("pheatmap")


#���ð�
library(edgeR)
library(pheatmap)

logFCfiler=2                  #logFC��������
fdrFilter=0.05                #�������pֵ��������
conNum=353                    #Ұ������Ʒ��Ŀ
treatNum=22                   #ͻ������Ʒ��Ŀ

setwd("C:\\Users\\lexb4\\Desktop\\snpRNA\\11.diff")                    #���ù���Ŀ¼
rt=read.table("sampleExp.txt",sep="\t",header=T,check.names=F)         #��ȡ�����ļ�

#���һ��������ֶ��У�ȡ��ֵ
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>1,]

#�����ݽ��н������������
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

#���������
write.table(diff,file="all.xls",sep="\t",quote=F)
diffSig = diff[(diff$FDR < fdrFilter & (diff$logFC>logFCfiler | diff$logFC<(-logFCfiler))),]
diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut, file="diff.xls",sep="\t",quote=F,col.names=F)
write.table(diffSigOut, file="diff.txt",sep="\t",quote=F,col.names=F)

#�����������ļ�
normalizeExp=rbind(id=colnames(newData),newData)
write.table(normalizeExp,file="normalExp.txt",sep="\t",quote=F,col.names=F)   #������л���У����ı���ֵ
diffExp=rbind(id=colnames(newData),newData[rownames(diffSig),])
write.table(diffExp,file="diffExp.txt",sep="\t",quote=F,col.names=F)         #����������У����ı���ֵ

#���Ʋ��������ͼ
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

#���ƻ�ɽͼ
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

