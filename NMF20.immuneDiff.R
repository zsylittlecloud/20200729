######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggpubr")


#引用包
library(limma)
library(ggpubr)

cluFile="cluster.txt"                #分型的结果文件
immFile="MCPcounter.result.txt"      #免疫细胞浸润结果文件
setwd("C:/Users/zsy/Desktop/巨噬细胞单细胞测序GBM/35MCP")     #设置工作目录

#读取免疫细胞浸润结果文件，并对数据进行整理
immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)
immune=as.matrix(immune)
immune=t(immune)

#读取分型的结果文件
cluster=read.table(cluFile, header=T, sep="\t", check.names=F, row.names=1)

#样品取交集
sameSample=intersect(row.names(immune), row.names(cluster))
immune1=immune[sameSample,,drop=F]
cluster1=cluster[sameSample,,drop=F]
data=cbind(immune1, cluster1)

#设置比较组
type=levels(factor(data[,"Cluster"]))
data$Cluster=factor(data$Cluster, levels=type)
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#设置颜色
bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:2]

#对免疫细胞进行循环
for(i in colnames(data)[1:(ncol(data)-1)]){
	#绘制小提琴图
	data[,i][data[,i]>quantile(data[,i],0.99)]=quantile(data[,i],0.99)
	violin=ggviolin(data, x="Cluster", y=i, fill = "Cluster", 
	         xlab="", ylab=i,
	         legend.title="Cluster",
	         palette=bioCol,
	         add = "boxplot", add.params = list(fill="white"))+ 
	         stat_compare_means(comparisons = my_comparisons)
	         #stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
	#输出图形
	pdf(file=paste0("violin.", i, ".pdf"), width=5, height=4.5)
	print(violin)
	dev.off()
}


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

