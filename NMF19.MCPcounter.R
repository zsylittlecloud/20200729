#引用包
library(limma)
library(MCPcounter)

expFile="symbol.txt"     #表达数据文件
setwd("C:/Users/zsy/Desktop/巨噬细胞单细胞测序GBM/35MCP")    #设置工作目录

#读取表达输入文件，并对输入文件处理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
mat=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
mat=avereps(mat)
data=mat[rowMeans(mat)>0,]

#删掉正常样品
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]
data=t(data)
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
data=t(avereps(data))

#免疫细胞浸润
MCPcounter.estimate=MCPcounter.estimate(data,
	featuresType="HUGO_symbols",
	probesets=read.table("probesets.txt",sep="\t",stringsAsFactors=FALSE,colClasses="character"),
	genes=read.table("genesets.txt",sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses="character",check.names=FALSE)
)
out=rbind(ID=colnames(MCPcounter.estimate), MCPcounter.estimate)
write.table(out, file="MCPcounter.result.txt", sep="\t", quote=F, col.names=F)


