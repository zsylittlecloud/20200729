#引用包
library(limma)
library(edgeR)
library(pheatmap)
library(ggplot2)

expFile="symbol.txt"      #表达输入文件
logFCfilter=1             #logFC过滤条件 
fdrFilter=0.05            #矫正后的p值过滤条件 
#pValue=0.05           #矫正前的p值过滤条件 
conNum=1156                 #前几个是正常样品数目 胶质瘤 5 169
treatNum=166             #后几个是肿瘤样品数目
setwd("C:/Users/zsy/Desktop/巨噬细胞单细胞测序GBM/1差异分析")      #设置工作目录

#读取数据文件，并对输入文件整理
rt=read.table(expFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#对数据进行矫正，数据过滤
dge=DGEList(counts=data)
cpmValue=cpm(dge)    #CPM是将counts转变为CPM指数counts
keep=rowMeans(cpmValue)>1
y=dge[keep, , keep.lib.sizes=FALSE]
y=calcNormFactors(y)
logCPM=cpm(y, log=TRUE, prior.count=3)  #logCPM值是根据CPM值通过log2(CPM + 2/L)计算得到的

#差异分析
Type=c(rep("con",conNum),rep("treat",treatNum))
design <- model.matrix(~0+factor(Type))
colnames(design) <- c("con","treat")
fit <- lmFit(logCPM,design)
cont.matrix<-makeContrasts(treat-con,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2,trend=TRUE)

#输出所有基因的差异结果
allDiff=topTable(fit2,adjust='fdr',number=200000)
write.table(allDiff,file="all.xls",sep="\t",quote=F)

#输出差异结果
diffSig <- allDiff[with(allDiff, (abs(logFC)>logFCfilter & adj.P.Val < fdrFilter )), ]   #P.Value<pValue
diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut,file="diff.xls",sep="\t",quote=F,col.names=F)
write.table(diffSigOut,file="diff.txt",sep="\t",quote=F,col.names=F)

#输出差异基因表达
heatmap=rbind(ID=colnames(data[as.vector(diffSigOut[,1]),]),data[as.vector(diffSigOut[,1]),])
write.table(heatmap,file="diffExp.txt",sep="\t",col.names=F,quote=F)

#绘制差异基因热图
geneNum=50
diffSig=diffSig[order(as.numeric(as.vector(diffSig$logFC))),]
diffGeneName=as.vector(rownames(diffSig))
diffLength=length(diffGeneName)
hmGene=c()
if(diffLength>(geneNum*2)){
  hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
}else{
  hmGene=diffGeneName
}
hmExp=logCPM[hmGene,]
Type=c(rep("N",conNum),rep("T",treatNum))
names(Type)=colnames(logCPM)
Type=as.data.frame(Type)
pdf(file="TCGA_heatmap.pdf",height=8,width=10)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("green", "black", "red"))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 8,
         fontsize_row=6,
         fontsize_col=8)
dev.off()


#定义显著性
Significant=ifelse((allDiff$adj.P.Val<fdrFilter & abs(allDiff$logFC)>logFCfilter), ifelse(allDiff$logFC>logFCfilter,"Up","Down"), "Not")
#绘制火山图
p = ggplot(allDiff, aes(logFC, -log10(adj.P.Val)))+
  geom_point(aes(col=Significant))+
  scale_color_manual(values=c("green", "black", "red"))+
  labs(title = " ")+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
p=p+theme_bw()
#保存为图片
pdf("TCGA_vol.pdf",width=5.5,height=5)
print(p)
dev.off()


#输出RPKM值 需geneLength.txt  RPKM更接近GEO数据可合并
#geneRT=read.table("geneLength.txt",sep="\t",header=T,check.names=F)
#geneRT=geneRT[,c(1,2,2)]
#geneRT=as.matrix(geneRT)
#rownames(geneRT)=geneRT[,1]
#exp1=geneRT[,2:ncol(geneRT)]
#dimnames1=list(rownames(exp1),colnames(exp1))
#data1=matrix(as.numeric(as.matrix(exp1)),nrow=nrow(exp1),dimnames=dimnames1)
#data1=avereps(data1)
#same <- intersect(row.names(y$counts),rownames(data1))
#data1=data1[same,]
#geneLen=data1[,1]
#rpkm=rpkm(y,geneLen)
#rpkm=rbind(id=colnames(rpkm),rpkm)
#write.table(rpkm,file="TCGA_rpkm.txt",sep="\t",quote=F,col.names=F)