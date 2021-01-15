#install.packages("rms")

library(rms)
#setwd("D:\\biowolf\\metabolism\\21.Nomogram")                         #???ù???Ŀ¼

#TCGA????ͼ????
riskFile="risk.txt"
cliFile="clinical-3.txt"
outFile="tcga.Nomogram.pdf"
risk=read.table(riskFile,header=T,sep="\t",check.names=F,row.names=1)        #??ȡ?????ļ?
cli=read.table(cliFile,sep="\t",check.names=F,header=T,row.names=1)          #??ȡ?ٴ??ļ?
sameSample=intersect(row.names(cli),row.names(risk))
risk=risk[sameSample,]
cli=cli[sameSample,]
rt=cbind(futime=risk[,1],fustat=risk[,2],cli,riskScore=risk[,(ncol(risk)-1)])
#???ݴ???
dd <- datadist(rt)
options(datadist="dd")
#???ɺ???
f <- cph(Surv(futime, fustat) ~ ., x=T, y=T, surv=T, data=rt, time.inc=1)
surv <- Survival(f)
#??��nomogram
nom <- nomogram(f, fun=list(function(x) surv(1, x), function(x) surv(2, x), function(x) surv(3, x)), 
    lp=F, funlabel=c("1-year survival", "2-year survival", "3-year survival"), 
    maxscale=100, 
    fun.at=c(0.99, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3,0.2,0.1,0.05))  
#nomogram???ӻ?
pdf(file=outFile,height=7.5,width=11)
plot(nom)
dev.off()

