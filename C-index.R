#???ð?
library(dplyr)
library(survival)
library(rms)
library(pec)

riskFile="risk.TCGA.txt"     #?????ļ?
cliFile="clinical.txt"      #?ٴ??????ļ?
setwd("C:/Users/zsy/Desktop/????ϸ????ϸ??????GBM/29C-index")     #???ù???Ŀ¼

#??ȡ?????ļ?
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk=risk[,c("futime", "fustat", "riskScore")]

#??ȡ?ٴ??????ļ?
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

#?ϲ?????
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1, cli)

#??????ɫ
bioCol=rainbow(ncol(rt)-1, s=0.9, v=0.9)

#????C-indexֵ
riskScore=cph(Surv(futime,fustat)~riskScore, data=rt, surv=TRUE)
Age=cph(Surv(futime,fustat)~Age, data=rt, surv=TRUE)
Gender=cph(Surv(futime,fustat)~Gender, data=rt, surv=TRUE)
Race=cph(Surv(futime,fustat)~Race, data=rt, surv=TRUE)
IDH=cph(Surv(futime,fustat)~IDH, data=rt, surv=TRUE)
c_index  <- cindex(list("Risk score"=riskScore, 
                        "Age"=Age,
                        "Gender"=Gender,
                        "Race"=Race,
                        "IDH"=IDH),
                    formula=Surv(futime,fustat)~ .,
                    data=rt,
                    eval.times=seq(0,10,1),
                    splitMethod="bootcv",
                    B=1000
                    )
#????ͼ??
pdf(file="C-index.pdf", width=5.5, height=5)
plot(c_index, 
     xlim=c(0,10), ylim=c(0.4,0.8), 
     col=bioCol, xlab="Time (years)",
     legend.x=6, legend.y=0.82, legend.cex=1)
dev.off()



