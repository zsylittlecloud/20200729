library(limma)                #���ð�
expFile="CGGA.share.txt"       #���������ļ�
cliFile="time.txt"            #���������ļ�
setwd("C:/Users/zsy/Desktop/����ϸ����ϸ������GBM/7CGGA�ϲ������ļ�")    #���ù���Ŀ¼

#��ȡ�����ļ������������ļ�����
rt=read.table(expFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
data=t(data)

#��ȡ��������
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)     #��ȡ�ٴ��ļ�

#���ݺϲ�
sameSample=intersect(row.names(data),row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]
out=cbind(cli,data)

#����ϲ���Ľ��
out=cbind(id=row.names(out),out)
write.table(out,file="CGGA.expTime.txt",sep="\t",row.names=F,quote=F)
