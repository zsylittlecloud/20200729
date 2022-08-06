######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("ggplot2")
#install.packages("ggalluvial")


#引用包
library(ggalluvial)
library(ggplot2)
library(dplyr)

clusterFile="cluster.txt"         #分型结果文件
subtypeFile="Subtype_Immune_Model_Based.txt"       #免疫分型文件
setwd("C:/Users/zsy/Desktop/巨噬细胞单细胞测序GBM/36ggalluvial")     #设置工作目录

#读取分型结果文件
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)

#读取免疫分型文件
subtype=read.table(subtypeFile, header=T, sep="\t", check.names=F, row.names=1)
rownames(subtype)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(subtype))

#样品取交集
sameSample=intersect(row.names(subtype), row.names(cluster))
subtype=subtype[sameSample,,drop=F]
subtype[,1]=gsub(".+\\((.+?)\\)", "\\1", subtype[,1])
colnames(subtype)=c("Immune subtype")
cluster=cluster[sameSample,,drop=F]
data=cbind(cluster, subtype)

#准备桑基图输入文件
corLodes=to_lodes_form(data, axes=1:ncol(data), id = "Cohort")

#得到输出文件
pdf(file="ggalluvial.pdf", width=5, height=7)
mycol=rep(c("#0066FF","#FF9900","#FF0000","#029149","#6E568C","#E0367A","#D8D155","#223D6C","#D20A13","#431A3D","#91612D","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#64495D","#7CC767"),15)
ggplot(corLodes, aes(x = x, stratum = stratum, alluvium = Cohort,fill = stratum, label = stratum)) +
  	 scale_x_discrete(expand = c(0, 0)) +  
  	 #用aes.flow控制线条颜色，forward说明颜色和前面的柱状图一致，backward说明和后面的柱状图一致。
  	 geom_flow(width = 2/10,aes.flow = "backward") + 
	 geom_stratum(alpha = .9,width = 2/10) +
	 scale_fill_manual(values = mycol) +
	 #size=3代表字体大小
	 geom_text(stat = "stratum", size = 3,color="black") +
	 xlab("") + ylab("") + theme_bw() + 
	 theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text.y = element_blank()) + #去掉坐标轴
	 theme(panel.grid =element_blank()) + 
	 theme(panel.border = element_blank()) + 
	 ggtitle("") + guides(fill = FALSE)                            
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

