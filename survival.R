#引用包
library(survival)
library(survminer)
setwd("C:/Users/zsy/Desktop/巨噬细胞单细胞测序GBM/10高低风险组生存分析")       #设置工作目录

#定义生存分析的函数
bioSurvival=function(inputFile=null, outFile=null){
	#读取输入文件
	rt=read.table(inputFile, header=T, sep="\t", check.names=F)
	
	#比较高低风险组生存差异，得到显著性p值
	diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
	pValue=1-pchisq(diff$chisq, df=1)
	if(pValue<0.001){
		pValue="p<0.001"
	}else{
		pValue=paste0("p=",sprintf("%.03f",pValue))
	}
	fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
		
	#绘制生存曲线
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

#调用函数，绘制生存曲线
bioSurvival(inputFile="risk.TCGA.txt", outFile="survival.TCGA.pdf")
bioSurvival(inputFile="risk.CGGA.txt", outFile="survival.CGGA.pdf")



