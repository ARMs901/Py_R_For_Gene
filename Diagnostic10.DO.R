#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")
#BiocManager::install("GSEABase")
#引用包
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library("GSEABase")
library("DOSE")

pvalueFilter=0.05       #p值过滤条件
qvalueFilter=0.05       #矫正后的p值过滤条件
#定义颜色
colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}
setwd("D:\\BaiduNetdiskDownload\\Diagnostic_video\\136Diagnostic\\10.DO")            #设置工作目录
rt=read.table("diff.txt", header=T, sep="\t", check.names=F)     #读取输入文件

#基因名字转换为基因id
genes=as.vector(rt[,1])
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=cbind(rt,entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]     #去除基因id为NA的基因

#DO富集分析
kk=enrichDO(gene=gene, ont="DO", pvalueCutoff=1, qvalueCutoff=1, readable=TRUE)
#使用enrichDO函数进行Disease Ontology (DO) 富集分析。参数中，指定用于分析的基因ID列表，
#设置p值和q值的初始过滤条件都为1（即不进行过滤），并要求结果是可读的。
DO=as.data.frame(kk)				#将富集分析的结果转换为数据框格式。
DO=DO[(DO$pvalue<pvalueFilter & DO$qvalue<qvalueFilter),]
									#根据设定的p值和q值过滤条件筛选出显著的DO富集结果。
#保存富集结果
write.table(DO, file="DO.txt", sep="\t", quote=F, row.names = F)

#定义显示疾病数目
showNum=30
if(nrow(DO)<showNum){
  showNum=nrow(DO)
}
#柱状图
pdf(file="barplot.pdf", width=6, height=11)
barplot(kk, drop=TRUE, showCategory=showNum, color=colorSel)
dev.off()
#气泡图
pdf(file="bubble.pdf", width = 6, height = 11)
dotplot(kk, showCategory=showNum, orderBy="GeneRatio", color=colorSel)
dev.off()
