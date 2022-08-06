#引用包
set.seed(123456)
library("glmnet")
library("survival")

coxSigFile="uniSigExp.txt"      #TCGA单因素显著基因的表达文件
geoFile="CGGA.expTime.txt"       #GEO数据库输入文件
setwd("C:/Users/zsy/Desktop/巨噬细胞单细胞测序GBM/9lasso和多因素cox模型建立")      #设置工作目录

#读取文件，并对输入文件进行整理
rt=read.table(coxSigFile, header=T, sep="\t", check.names=F, row.names=1)
rt$futime[rt$futime<=0]=0.003

#构建lasso回归模型
x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$futime, rt$fustat))
fit=glmnet(x, y, family="cox", maxit=1000)
#lasso回归图形
pdf("lasso.lambda.pdf")
plot(fit, xvar="lambda", label=TRUE)
dev.off()
#交叉验证图形
cvfit=cv.glmnet(x, y, family="cox", maxit=1000)
pdf("lasso.cvfit.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)), lty="dashed")
dev.off()

#输出lasso显著基因表达量
coef=coef(fit, s=cvfit$lambda.min)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
lassoSigExp=rt[,c("futime", "fustat", lassoGene)]
lassoSigExpOut=cbind(id=row.names(lassoSigExp), lassoSigExp)
write.table(lassoSigExpOut,file="lasso.SigExp.txt",sep="\t",row.names=F,quote=F)

#构建COX模型
multiCox=coxph(Surv(futime, fustat) ~ ., data = lassoSigExp)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)

#输出模型的公式
outMultiTab=data.frame()
outMultiTab=cbind(
		          coef=multiCoxSum$coefficients[,"coef"],
		          HR=multiCoxSum$conf.int[,"exp(coef)"],
		          HR.95L=multiCoxSum$conf.int[,"lower .95"],
		          HR.95H=multiCoxSum$conf.int[,"upper .95"],
		          pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outMultiTab=cbind(id=row.names(outMultiTab), outMultiTab)
outMultiTab=outMultiTab[,1:2]
write.table(outMultiTab, file="multiCox.txt", sep="\t", row.names=F, quote=F)

#输出train组风险值
trainScore=predict(multiCox, type="risk", newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
risk=as.vector(ifelse(trainScore>median(trainScore),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(trainScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="risk.TCGA.txt",sep="\t",quote=F,row.names=F)

#输出test组风险值
rt=read.table(geoFile, header=T, sep="\t", check.names=F, row.names=1)
rt$futime=rt$futime/365
testFinalGeneExp=rt[,coxGene]
testScore=predict(multiCox,type="risk",newdata=rt)
outCol=c("futime","fustat",coxGene)
risk=as.vector(ifelse(testScore>median(trainScore),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(testScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="risk.CGGA.txt",sep="\t",quote=F,row.names=F)



