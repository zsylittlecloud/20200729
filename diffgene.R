library(limma)
#设置工作目录
setwd("C:/Users/zsy/Desktop/巨噬细胞单细胞测序GBM/2提取差异的基因")
rt=read.table("symbol.txt",sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

gene=read.table("gene.txt", header=F, check.names=F, sep="\t")
sameGene=intersect(as.vector(gene[,1]),rownames(data))
geneExp=data[sameGene,]

out=rbind(ID=colnames(geneExp),geneExp)
write.table(out,file="symbolgene.txt",sep="\t",quote=F,col.names=F)