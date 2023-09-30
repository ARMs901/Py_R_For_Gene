# ��������Ŀ�
library(limma)
library(ggpubr)

# ���������ļ���·��
expFile="U_D_GSE24287.txt"          
conFile="U_D_GSE24287_s1.txt"      
treatFile="U_D_GSE24287_s2.txt"     
geneFile="intersected_genes.txt"     
setwd("D:\\BaiduNetdiskDownload\\Diagnostic_video\\136Diagnostic\\15.testDiff\\U_D") # ���ù���Ŀ¼

# ��ȡ���������ļ�
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]

# ���¶������ݵ�ά��
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
rt=avereps(data)

# �Ը�ֵ����ȡ����
qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
if(LogC){
	rt[rt<0]=0
	rt=log2(rt+1)}

# ��׼������
data=normalizeBetweenArrays(rt)

# ��ȡ������Ϣ������
con=read.table(conFile, header=F, sep="\t", check.names=F)
treat=read.table(treatFile, header=F, sep="\t", check.names=F)
conData=data[,as.vector(con[,1])]
treatData=data[,as.vector(treat[,1])]
data=cbind(conData, treatData)
conNum=ncol(conData)
treatNum=ncol(treatData)

# �����׼���������
Type=c(rep("con",conNum), rep("treat",treatNum))
outData=rbind(id=paste0(colnames(data),"_",Type),data)
write.table(outData, file="test.normalize.txt", sep="\t", quote=F, col.names=F)

# �ӽ��������ļ��ж�ȡ����
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),,drop=F]

# ���ñȽ���
Type=c(rep("Con",conNum), rep("Treat",treatNum))
my_comparisons=list()
my_comparisons[[1]]=levels(factor(Type))

# ��ÿ��������в���������������������ͼ
newGeneLists=c()
outTab=data.frame()
for(i in row.names(data)){
	# �Ա������ݽ��е����������쳣�ߵ�ֵ
	rt1=data.frame(expression=data[i,], Type=Type)
	
	# ��������ͼ�����в���Ƚ�
	boxplot=ggboxplot(rt1, x="Type", y="expression", color="Type",
				      xlab="",
				      ylab=paste(i, "expression"),
				      legend.title="",
				      palette = c("blue", "red"),
				      add = "jitter")+ 
		stat_compare_means(comparisons = my_comparisons)
		
	# �������ͼ��PDF�ļ�
	pdf(file=paste0("boxplot.",i,".pdf"), width=5, height=4.5)
	print(boxplot)
	dev.off()
}