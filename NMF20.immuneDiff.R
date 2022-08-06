######Video source: https://ke.biowolf.cn
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺biowolf@foxmail.com
######����΢��: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggpubr")


#���ð�
library(limma)
library(ggpubr)

cluFile="cluster.txt"                #���͵Ľ���ļ�
immFile="MCPcounter.result.txt"      #����ϸ���������ļ�
setwd("C:/Users/zsy/Desktop/����ϸ����ϸ������GBM/35MCP")     #���ù���Ŀ¼

#��ȡ����ϸ���������ļ����������ݽ�������
immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)
immune=as.matrix(immune)
immune=t(immune)

#��ȡ���͵Ľ���ļ�
cluster=read.table(cluFile, header=T, sep="\t", check.names=F, row.names=1)

#��Ʒȡ����
sameSample=intersect(row.names(immune), row.names(cluster))
immune1=immune[sameSample,,drop=F]
cluster1=cluster[sameSample,,drop=F]
data=cbind(immune1, cluster1)

#���ñȽ���
type=levels(factor(data[,"Cluster"]))
data$Cluster=factor(data$Cluster, levels=type)
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#������ɫ
bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:2]

#������ϸ������ѭ��
for(i in colnames(data)[1:(ncol(data)-1)]){
	#����С����ͼ
	data[,i][data[,i]>quantile(data[,i],0.99)]=quantile(data[,i],0.99)
	violin=ggviolin(data, x="Cluster", y=i, fill = "Cluster", 
	         xlab="", ylab=i,
	         legend.title="Cluster",
	         palette=bioCol,
	         add = "boxplot", add.params = list(fill="white"))+ 
	         stat_compare_means(comparisons = my_comparisons)
	         #stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
	#���ͼ��
	pdf(file=paste0("violin.", i, ".pdf"), width=5, height=4.5)
	print(violin)
	dev.off()
}


######Video source: https://ke.biowolf.cn
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺biowolf@foxmail.com
######����΢��: 18520221056
