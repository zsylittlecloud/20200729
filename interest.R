library(venn)                  
outFile="intersectGenes.txt"                      #输出文件
setwd("C:/Users/zsy/Desktop/巨噬细胞单细胞测序GBM/4交集韦恩图")    
geneList=list()

rt=read.table("Markgene.txt",sep="\t",header=F,check.names=F)     #输入文件第一列为基因且有表头
geneNames=as.vector(rt[,1])              
geneNames=gsub("^ | $","",geneNames)     #对于gsub函数，简单理解为gsub(“文本中要处理的特殊符号”，”替换为的符号“，文件名)  
uniqGene=unique(geneNames)                   #去重           
geneList[["Markgene"]]=uniqGene                      #类别命名

rt=read.table("TCGA_turquoise brown.txt",sep="\t",header=F,check.names=F)
geneNames=as.vector(rt[,1])               
geneNames=gsub("^ | $","",geneNames)    
uniqGene=unique(geneNames)             
geneList[["TCGA_turquoise brown"]]=uniqGene


mycol=c("#029149","#E0367A","#5D90BA","#431A3D","#FFD121","#D8D155","#223D6C","#D20A13","#088247","#11AA4D","#7A142C","#5D90BA","#64495D","#7CC767")
pdf(file="venn.pdf",width=5,height=5)
venn(geneList,col=mycol[1:length(geneList)],zcolor=mycol[1:length(geneList)],box=F,ilabels=F)
dev.off()

intersectGenes=Reduce(intersect,geneList)
write.table(file=outFile,intersectGenes,sep="\t",quote=F,col.names=F,row.names=F)
