#引用包
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

pvalueFilter=0.05       #p值过滤条件
qvalueFilter=0.05       #矫正后的p值过滤条件

#定义颜色
colorSel="qvalue"		#默认颜色选择基于q值。
if(qvalueFilter>0.05){
  colorSel="pvalue"
}						#如果q值过滤条件大于0.05，则颜色选择基于p值。

setwd("D:\\BaiduNetdiskDownload\\Diagnostic_video\\136Diagnostic\\08.GO")           #设置工作目录
rt=read.table("diff.txt", header=T, sep="\t", check.names=F)     #读取输入文件

#基因名字转换为基因id
genes=as.vector(rt[,1])					#获取第一列的数据作为基因名。
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)		#将结果转换为字符向量。
gene=entrezIDs[entrezIDs!="NA"]        	#移除任何为NA的Entrez IDs。

#GO富集分析
kk=enrichGO(gene=gene, OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable=T)
GO=as.data.frame(kk)					#将结果转换为数据框。
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
										#保存富集结果，仅保留满足特定p值和q值条件的结果。
write.table(GO, file="GO.txt", sep="\t", quote=F, row.names = F)
										#将过滤后的富集结果保存到"GO.txt"文件。
#定义显示GO的数目
showNum=10								#默认显示10个类别。
if(nrow(GO)<30){
  showNum=nrow(GO)
}										#如果结果中的GO类别少于30个，则全部显示。


#柱状图
pdf(file="barplot.pdf", width=7*0.9*2, height=11*0.8*2)	#设定输出文件为PDF
bar=barplot(kk, drop=TRUE, showCategory=showNum, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()

#气泡图
pdf(file="bubble.pdf", width=7*0.9*2, height=11*0.8*2)
bub=dotplot(kk, showCategory=showNum, orderBy="GeneRatio", split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bub)
dev.off()

