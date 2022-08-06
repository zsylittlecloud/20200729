#���ð�
library(limma)
library(reshape2)
library(ggpubr)

riskFile="risk.TCGA.txt"       #�����ļ�
scoreFile="TMEscores.txt"      #����΢��������ļ�
setwd("C:/Users/zsy/Desktop/����ϸ����ϸ������GBM/22��ֺ�TME")     #���ù���Ŀ¼

#��ȡ�����ļ�
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk$risk=factor(risk$risk, levels=c("low","high"))

#��ȡ����΢��������ļ����������ݽ�������
score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
score=as.matrix(score)
row.names(score)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(score))
score=avereps(score)
score=score[,1:3]

#��Ʒȡ����
sameSample=intersect(row.names(risk), row.names(score))
risk=risk[sameSample,"risk",drop=F]
score=score[sameSample,,drop=F]
rt=cbind(risk, score)

#���ϲ��������ת��Ϊggplot2�������ļ�
data=melt(rt, id.vars=c("risk"))
colnames(data)=c("Risk", "scoreType", "Score")

#����С����ͼ
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

#���ͼ��
pdf(file="vioplot.pdf", width=6, height=5)
print(p1)
dev.off()

