#���ð�
library(limma)
library(MCPcounter)

expFile="symbol.txt"     #���������ļ�
setwd("C:/Users/zsy/Desktop/����ϸ����ϸ������GBM/35MCP")    #���ù���Ŀ¼

#��ȡ���������ļ������������ļ�����
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
mat=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
mat=avereps(mat)
data=mat[rowMeans(mat)>0,]

#ɾ��������Ʒ
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]
data=t(data)
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
data=t(avereps(data))

#����ϸ������
MCPcounter.estimate=MCPcounter.estimate(data,
	featuresType="HUGO_symbols",
	probesets=read.table("probesets.txt",sep="\t",stringsAsFactors=FALSE,colClasses="character"),
	genes=read.table("genesets.txt",sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses="character",check.names=FALSE)
)
out=rbind(ID=colnames(MCPcounter.estimate), MCPcounter.estimate)
write.table(out, file="MCPcounter.result.txt", sep="\t", quote=F, col.names=F)

