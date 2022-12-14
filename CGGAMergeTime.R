library(limma)                #引用包
expFile="CGGA.share.txt"       #表达数据文件
cliFile="time.txt"            #生存数据文件
setwd("C:/Users/zsy/Desktop/巨噬细胞单细胞测序GBM/7CGGA合并生存文件")    #设置工作目录

#读取表达文件，并对输入文件整理
rt=read.table(expFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
data=t(data)

#读取生存数据
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)     #读取临床文件

#数据合并
sameSample=intersect(row.names(data),row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]
out=cbind(cli,data)

#输出合并后的结果
out=cbind(id=row.names(out),out)
write.table(out,file="CGGA.expTime.txt",sep="\t",row.names=F,quote=F)

