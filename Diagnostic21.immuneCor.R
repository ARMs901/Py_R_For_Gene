
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("reshape2")
#install.packages("ggpubr")
#install.packages("ggExtra")


#引用包
library(limma)
library(reshape2)
library(ggpubr)
library(ggExtra)

gene="COL6A3"                       #基因名称
expFile="normalize.txt"             #表达数据文件
immFile="CIBERSORT-Results.txt"     #免疫细胞浸润结果文件
setwd("D:\\BaiduNetdiskDownload\\Diagnostic_video\\136Diagnostic\\21.immuneCor\\U_D")     #设置工作目录

#读取基因表达文件,并对数据进行处理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#获取目标基因表达量
data=t(data[gene,,drop=F])
data=as.data.frame(data)

#读取免疫细胞结果文件，并对数据进行整理
immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)

#将基因表达数据和免疫细胞浸润数据合并
sameSample=intersect(row.names(immune), row.names(data))
rt=cbind(immune[sameSample,,drop=F], data[sameSample,,drop=F])

# 初始化输出数据框
outTab=data.frame()

# 对每种免疫细胞，计算其与目标基因的相关性，并绘制散点图
for(i in colnames(rt)[1:(ncol(rt)-1)]){
	x=as.numeric(rt[,gene])
	y=as.numeric(rt[,i])
	if(sd(y)==0){y[1]=0.00001}
	cor=cor.test(x, y, method="spearma")
	
	outVector=cbind(Gene=gene, Cell=i, cor=cor$estimate, pvalue=cor$p.value)
	outTab=rbind(outTab,outVector)
	
	if(cor$p.value<0.05){
		outFile=paste0("cor.", i, ".pdf")
		df1=as.data.frame(cbind(x,y))
		p1=ggplot(df1, aes(x, y)) + 
				  xlab(paste0(gene, " expression")) + ylab(i)+
				  geom_point() + geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
				  stat_cor(method = 'spearman', aes(x =x, y =y))
		p2=ggMarginal(p1, type="density", xparams=list(fill = "orange"), yparams=list(fill = "blue"))
		#相关性图形
		pdf(file=outFile, width=5.2, height=5)
		print(p2)
		dev.off()
	}
}
#输出免疫功能和p值表格文件
write.table(outTab,file="cor.result.txt",sep="\t",row.names=F,quote=F)


