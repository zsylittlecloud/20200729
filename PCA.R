#���ð�
library(limma)
library(scatterplot3d)
setwd("C:/Users/zsy/Desktop/����ϸ����ϸ������GBM/28PCA")        #���ù���Ŀ¼

#����PCA��������
myPCA=function(input=null,output=null){
	#��ȡ���������ļ�
	rt=read.table(input, header=T, sep="\t", check.names=F)
	rt=as.matrix(rt)
	rownames(rt)=rt[,1]
	exp=rt[,2:ncol(rt)]
	dimnames=list(rownames(exp),colnames(exp))
	data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
	data=avereps(data)
	data=data[rowMeans(data)>0.5,]
	
	#ɾ��������Ʒ
	#type=sapply(strsplit(colnames(data),"\\-"),"[",4)
	#type=sapply(strsplit(type,""),"[",1)
	#type=gsub("2","1",type)
	data=t(data)
	#rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",rownames(data))
		
	#��ȡrisk�����ļ�
	risk=read.table("risk.CGGA.txt", header=T, sep="\t", row.names=1, check.names=F)
	sameSample=intersect(rownames(data),rownames(risk))
	data=data[sameSample,]
	risk=risk[sameSample,]
	group=as.vector(risk[,"risk"])
		
	#PCA����
	data.class <- rownames(data)
	data.pca <- prcomp(data, scale. = TRUE)

	#����PCAͼ��
	color=ifelse(group=="low",4,2)
	pcaPredict=predict(data.pca)
	pdf(file=output, width=7, height=7)
	par(oma=c(1,1,2.5,1))
	s3d=scatterplot3d(pcaPredict[,1:3], pch = 16, color=color, angle=35)
	legend("top", legend = c("Low risk","High risk"),pch = 16, inset = -0.2, box.col="white", xpd = TRUE, horiz = TRUE,col=c(4,2))
	dev.off()
}

######�������л����PCAͼ����04�ڿ�symbol.txt���Ƶ���ǰĿ¼
myPCA(input="CGGAGBM.txt",output="PCA.allGene1.pdf")
######����ͭ������ػ����PCAͼ����09�ڿ�cuproptosisExp.txt���Ƶ���ǰĿ¼
#myPCA(input="cuproptosisExp.txt",output="PCA.cuproptosisGene.pdf")
######����ͭ�������lncRNA��PCAͼ����09�ڿ�cuproptosisLncExp.txt���Ƶ���ǰĿ¼
#myPCA(input="cuproptosisLncExp.txt",output="PCA.cuproptosisLncRNA.pdf")


######��ȡ�����ļ�,����ģ��lncRNA��PCAͼ����14�ڿ�risk.all.txt���Ƶ���ǰĿ¼
risk=read.table("risk.CGGA.txt", header=T, sep="\t", check.names=F, row.names=1)
data=risk[,3:(ncol(risk)-2)]
group=as.vector(risk[,"risk"])
		
#PCA����
data.class <- rownames(data)
data.pca <- prcomp(data, scale. = TRUE)

#���ӻ�
color=ifelse(group=="low",4,2)
pcaPredict=predict(data.pca)
pdf(file="PCA.riskLnc1.pdf", width=6.5, height=6)
par(oma=c(1,1,2.5,1))
s3d=scatterplot3d(pcaPredict[,1:3], pch = 16, color=color, angle=35)
legend("top", legend = c("Low risk","High risk"),pch = 16, inset = -0.2, box.col="white", xpd = TRUE, horiz = TRUE,col=c(4,2))
dev.off()

