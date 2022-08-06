library(limma)                #???ð?
expFile="tcga.share.txt"      #?????????ļ?
cliFile="time.txt"            #?????????ļ?
setwd("C:/Users/zsy/Desktop/巨噬细胞单细胞测序GBM/6TCGA合并生存文件")      #???ù???Ŀ¼

#??ȡ?????ļ????????????ļ?????
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
data=t(data)

#??ȡ????????
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

#???ݺϲ?
sameSample=intersect(row.names(data), row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]
out=cbind(cli, data)

#?????ϲ????Ľ???
out=cbind(id=row.names(out), out)
write.table(out, file="TCGA.expTime.txt", sep="\t", row.names=F, quote=F)


