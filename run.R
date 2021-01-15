#install.packages('e1071')

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("preprocessCore")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


library("limma")                                                      #���ð�
setwd("D:\\biowolf\\TMBimmune\\19.CIBERSORT")          #���ù���Ŀ¼
expFile="symbol.txt"                                                  #���������ļ�

#��ȡ�����ļ������������ļ�����
rt=read.table(expFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#ȥ��������Ʒ
if(grepl("-",colnames(data)[ncol(data)])){
	group=sapply(strsplit(colnames(data),"\\-"),"[",4)
	group=sapply(strsplit(group,""),"[",1)
	group=gsub("2","1",group)
	data=data[,group==0]
}

#ȥ���ͱ������
data=data[rowMeans(data)>0,]

#���ݽ���
v <-voom(data, plot = F, save.plot = F)
out=v$E
out=rbind(ID=colnames(out),out)
write.table(out,file="uniq.symbol.txt",sep="\t",quote=F,col.names=F)        #����ļ�

#����CIBERSORT���õ�����ϸ���������
source("TMBimmune19.CIBERSORT.R")
results=CIBERSORT("ref.txt", "uniq.symbol.txt", perm=100, QN=TRUE)

