#引用包
library(limma)
library(reshape2)
library(ggplot2)

expFile="TCGA.normalize.txt"    #表达数据文件
geneFile="gene.txt"             #免疫检查点基因列表文件
riskFile="risk.TCGA.txt"        #风险文件
setwd("C:/Users/zsy/Desktop/巨噬细胞单细胞测序GBM/23免疫检查点相关性")      #设置工作目录

#读取表达输入文件,并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

#读取基因列表文件，获取免疫检查点相关基因的表达量
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(geneRT[,1]), rownames(data))
data=t(data[sameGene,])

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(data), row.names(risk))
data=data[sameSample,,drop=F]
risk=risk[sameSample,3:(ncol(risk)-1),drop=F]

#相关性分析
outTab=data.frame()
for(checkpiont in colnames(data)){
	for(gene in colnames(risk)){
		x=as.numeric(data[,checkpiont])
		y=as.numeric(risk[,gene])
		corT=cor.test(x,y,method="spearman")
		cor=corT$estimate
		pvalue=corT$p.value
		text=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
		outTab=rbind(outTab,cbind(Gene=gene, checkpiont=checkpiont, cor, text, pvalue))
	}
}

#绘制相关性热图
outTab$Gene=factor(outTab$Gene, levels=colnames(risk))
outTab$cor=as.numeric(outTab$cor)
pdf(file="checkpointCor.pdf", width=8, height=7)
ggplot(outTab, aes(Gene, checkpiont)) + 
	geom_tile(aes(fill = cor), colour = "grey", size = 1)+
	scale_fill_gradient2(low = "#5C5DAF", mid = "white", high = "#EA2E2D") + 
	geom_text(aes(label=text),col ="black",size = 3) +
	theme_minimal() +    #去掉背景
	theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(),
	      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),   #x轴字体
	      axis.text.y = element_text(size = 10, face = "bold")) +       #y轴字体
	labs(fill =paste0("***  p<0.001","\n", "**  p<0.01","\n", " *  p<0.05","\n", "\n","Correlation")) +   #设置图例
	scale_x_discrete(position = "bottom")      #X轴名称显示位置
dev.off()


