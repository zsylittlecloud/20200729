library(ggpubr)              #引用包
tciaFile="TCIA.txt"          #TCIA打分文件
riskFile="risk.TCGA.txt"     #风险文件
setwd("C:/Users/zsy/Desktop/巨噬细胞单细胞测序GBM/27TCIA")     #修改工作目录

#读取TCIA打分文件
ips=read.table(tciaFile, header=T, sep="\t", check.names=F, row.names=1)

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
	
#合并数据
sameSample=intersect(row.names(ips), row.names(risk))
ips=ips[sameSample, , drop=F]
risk=risk[sameSample, "risk", drop=F]
data=cbind(ips, risk)
	
#设置比较组
data$risk=ifelse(data$risk=="high", "High-risk", "Low-risk")
group=levels(factor(data$risk))
data$risk=factor(data$risk, levels=c("Low-risk", "High-risk"))
group=levels(factor(data$risk))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#对TCIA打分进行循环,分别绘制小提琴图
for(i in colnames(data)[1:(ncol(data)-1)]){
	rt=data[,c(i, "risk")]
	colnames(rt)=c("IPS", "Risk")
	gg1=ggviolin(rt, x="Risk", y="IPS", fill = "Risk", 
	         xlab="", ylab=i,
	         legend.title="Risk",
	         add = "boxplot", add.params = list(fill="white"))+ 
	         stat_compare_means(comparisons = my_comparisons)
	         #stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
	
	pdf(file=paste0(i, ".pdf"), width=6, height=5)
	print(gg1)
	dev.off()
}



