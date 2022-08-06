library(ggpubr)              #���ð�
tciaFile="TCIA.txt"          #TCIA����ļ�
riskFile="risk.TCGA.txt"     #�����ļ�
setwd("C:/Users/zsy/Desktop/����ϸ����ϸ������GBM/27TCIA")     #�޸Ĺ���Ŀ¼

#��ȡTCIA����ļ�
ips=read.table(tciaFile, header=T, sep="\t", check.names=F, row.names=1)

#��ȡ�����ļ�
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
	
#�ϲ�����
sameSample=intersect(row.names(ips), row.names(risk))
ips=ips[sameSample, , drop=F]
risk=risk[sameSample, "risk", drop=F]
data=cbind(ips, risk)
	
#���ñȽ���
data$risk=ifelse(data$risk=="high", "High-risk", "Low-risk")
group=levels(factor(data$risk))
data$risk=factor(data$risk, levels=c("Low-risk", "High-risk"))
group=levels(factor(data$risk))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#��TCIA��ֽ���ѭ��,�ֱ����С����ͼ
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


