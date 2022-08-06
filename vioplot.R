#引用包
library(limma)
library(reshape2)
library(ggpubr)

riskFile="risk.TCGA.txt"       #风险文件
scoreFile="TMEscores.txt"      #肿瘤微环境打分文件
setwd("C:/Users/zsy/Desktop/巨噬细胞单细胞测序GBM/22打分和TME")     #设置工作目录

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk$risk=factor(risk$risk, levels=c("low","high"))

#读取肿瘤微环境打分文件，并对数据进行整理
score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
score=as.matrix(score)
row.names(score)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(score))
score=avereps(score)
score=score[,1:3]

#样品取交集
sameSample=intersect(row.names(risk), row.names(score))
risk=risk[sameSample,"risk",drop=F]
score=score[sameSample,,drop=F]
rt=cbind(risk, score)

#将合并后的数据转换为ggplot2的输入文件
data=melt(rt, id.vars=c("risk"))
colnames(data)=c("Risk", "scoreType", "Score")

#绘制小提琴图
p=ggviolin(data, x="scoreType", y="Score", fill = "Risk",
	     xlab="",
	     ylab="TME score",
	     legend.title="Risk",
	     add = "boxplot", add.params = list(color="white"),
	     palette = c("blue","red"), width=1)
p=p+rotate_x_text(45)
p1=p+stat_compare_means(aes(group=Risk),
	      method="wilcox.test",
	      symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
	      label = "p.signif")

#输出图形
pdf(file="vioplot.pdf", width=6, height=5)
print(p1)
dev.off()


