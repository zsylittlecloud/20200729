
#引用包
library(ggpubr)
library(reshape2)

riskFile="risk.TCGA.txt"      #风险文件
tmbFile="TMB.txt"             #肿瘤突变负荷文件
setwd("C:/Users/zsy/Desktop/巨噬细胞单细胞测序GBM/18肿瘤突变负荷相关性")      #设置工作目录

#读取输入文件
tmb=read.table(tmbFile, header=T, sep="\t", check.names=F, row.names=1)      #读取肿瘤突变负荷文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)    #读取风险文件

#合并数据
tmb=as.matrix(tmb)
tmb[tmb>quantile(tmb,0.975)]=quantile(tmb,0.975)
sameSample=intersect(row.names(tmb), row.names(risk))
tmb=tmb[sameSample,,drop=F]
risk=risk[sameSample,,drop=F]
data=cbind(risk, tmb)
data=data[,c("riskScore", "risk", "TMB")]

#设置比较组
data$risk=factor(data$risk, levels=c("low", "high"))
risk=levels(factor(data$risk))
comp=combn(risk, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#设置高低风险组的颜色
bioCol=c("#0066FF","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(risk)]

#绘制箱线图
boxplot=ggboxplot(data, x="risk", y="TMB", fill="risk",
		          xlab="",
		          ylab="Tumor Burden Mutation",
		          legend.title="Risk",
		          palette = bioCol )+ 
	    stat_compare_means(comparisons = my_comparisons)
pdf(file="boxplot.pdf",width=5,height=4.5)
print(boxplot)
dev.off()

#相关性图形
length=length(levels(factor(data$risk)))
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
p1=ggplot(data, aes(riskScore, TMB)) + 
		  xlab("Risk score")+ylab("Tumor Burden Mutation")+
		  geom_point(aes(colour=risk))+
		  scale_color_manual(values=bioCol[1:length])+ 
		  geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
		  stat_cor(method = 'spearman', aes(x =riskScore, y =TMB))
#相关性图形
pdf(file="cor.pdf", width=6, height=4.5)
print(p1)
dev.off()



