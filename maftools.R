library(maftools)       #引用包
setwd("C:/Users/zsy/Desktop/巨噬细胞单细胞测序GBM/17风险打分突变瀑布图")      #设置工作目录

#读取风险文件
risk=read.table("risk.TCGA.txt", header=T, sep="\t", check.names=F)
outTab=risk[,c(1, ncol(risk))]
colnames(outTab)=c("Tumor_Sample_Barcode", "Risk")
write.table(outTab, file="ann.txt", sep="\t", quote=F, row.names=F)

#读取基因突变文件
geneNum=20    #基因的数目，展示突变频率最高的20个基因
geneMut=read.table("geneMut.txt", header=T, sep="\t", check.names=F, row.names=1)
gene=row.names(geneMut)[1:geneNum]

#定义注释的颜色
ann_colors=list()
col=c("blue", "red")
names(col)=c("low", "high")
ann_colors[["Risk"]]=col

#低风险组瀑布图
pdf(file="low.pdf", width=6, height=6)
maf=read.maf(maf="low.maf", clinicalData="ann.txt")
oncoplot(maf=maf, clinicalFeatures="Risk", genes=gene, annotationColor=ann_colors, keepGeneOrder=T)
dev.off()

#高风险组瀑布图
pdf(file="high.pdf", width=6, height=6)
maf=read.maf(maf="high.maf", clinicalData="ann.txt")
oncoplot(maf=maf, clinicalFeatures="Risk", genes=gene, annotationColor=ann_colors, keepGeneOrder=T)
dev.off()



