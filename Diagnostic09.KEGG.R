
#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")
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

setwd("D:\\BaiduNetdiskDownload\\Diagnostic_video\\136Diagnostic\\09.KEGG")          #设置工作目录
rt=read.table("diff.txt", header=T, sep="\t", check.names=F)     #读取输入文件

#基因名字转换为基因id
genes=as.vector(rt[,1])	#从数据框rt的第一列获取基因名称，并转化为向量形式。
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
						#使用mget函数，尝试将基因名转换为对应的Entrez IDs。
						#如果找不到对应的ID，将返回NA
entrezIDs=as.character(entrezIDs)		#将获取到的Entrez IDs转换为字符向量
rt=cbind(rt,entrezID=entrezIDs)			#将Entrez IDs添加到原始数据框rt中，作为新的一列。
gene=entrezIDs[entrezIDs!="NA"]        	#去除基因id为NA的基因

#kegg富集分析
kk=enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
#使用enrichKEGG函数进行KEGG富集分析。这里的参数表示用于分析的基因ID列表，分析的生物体为人类（"hsa"），
#并设置p值和q值的初始过滤条件都为1（即不进行过滤）。
KEGG=as.data.frame(kk)
#将富集分析的结果转换为数据框格式。
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$id[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
#将KEGG结果中的基因ID转换回原始的基因名称，并将这些名称连接起来形成一个字符串。
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
#根据设定的p值和q值过滤条件筛选出显著的KEGG富集结果。
#保存富集结果
write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)

#定义显示通路数目
showNum=30
if(nrow(KEGG)<showNum){
  showNum=nrow(KEGG)
}

#柱状图
pdf(file="barplot.pdf", width=8, height=14)
barplot(kk, drop=TRUE, showCategory=showNum, color=colorSel)
dev.off()

#气泡图
pdf(file="bubble.pdf", width = 8, height = 14)
dotplot(kk, showCategory=showNum, orderBy="GeneRatio", color=colorSel)
dev.off()