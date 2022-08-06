#���ð�
library(survival)
library(survminer)
setwd("C:/Users/zsy/Desktop/����ϸ����ϸ������GBM/10�ߵͷ������������")       #���ù���Ŀ¼

#������������ĺ���
bioSurvival=function(inputFile=null, outFile=null){
	#��ȡ�����ļ�
	rt=read.table(inputFile, header=T, sep="\t", check.names=F)
	
	#�Ƚϸߵͷ�����������죬�õ�������pֵ
	diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
	pValue=1-pchisq(diff$chisq, df=1)
	if(pValue<0.001){
		pValue="p<0.001"
	}else{
		pValue=paste0("p=",sprintf("%.03f",pValue))
	}
	fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
		
	#������������
	surPlot=ggsurvplot(fit, 
		           data=rt,
		           conf.int=T,
		           pval=pValue,
		           pval.size=6,
		           legend.title="Risk",
		           legend.labs=c("High risk", "Low risk"),
		           xlab="Time(years)",
		           break.time.by = 1,
		           palette=c("red", "blue"),
		           risk.table=TRUE,
		           risk.table.title="",
		           risk.table.height=.25)
	pdf(file=outFile, width=6.5, height=5.5, onefile=FALSE)
	print(surPlot)
	dev.off()
}

#���ú�����������������
bioSurvival(inputFile="risk.TCGA.txt", outFile="survival.TCGA.pdf")
bioSurvival(inputFile="risk.CGGA.txt", outFile="survival.CGGA.pdf")


