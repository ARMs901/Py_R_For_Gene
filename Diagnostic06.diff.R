
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
#install.packages("pheatmap")

#���ð�
library(limma)
library(pheatmap)

inputFile="merge.txt"       #�����ļ�
logFCfilter=2               #logFC������ֵ
adj.P.Val.Filter=0.05       #������pֵ��ֵ
setwd("C:\\biowolf\\Diagnostic\\06.diff")      #���ù���Ŀ¼

#��ȡ�����ļ������������ļ�����
#ʹ��read.table������inputFile�ж�ȡ���ݡ����У�
#header=T���ļ��ĵ�һ����������
#sep="\t"���ֶ������Ʊ�����tab���ָ��ġ�
#check.names=F�����޸�������ʹ���Ϊ��׼��R�����������磬������������ո�������ַ���Rͨ�����޸����ǣ���
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
#��rtת��Ϊһ������
rt=as.matrix(rt)
#��rt����ĵ�һ������Ϊ������
rownames(rt)=rt[,1]
#��ȡrt����ĵڶ��е����һ�У��������������exp�С�
exp=rt[,2:ncol(rt)]
#Ϊ��һ�������ݾ��󴴽�ά�����ơ������ȡexp���������������������Ƿ���һ���б��У�������dimnames�С�
dimnames=list(rownames(exp),colnames(exp))
#��exp�ٴ�ת��Ϊһ����ֵ����ȷ������Ԫ�ض������֣�����Ϊ����¾�������������������
#��������Ϊ��ȷ�����е�Ԫ�ض�����ֵ���͡�
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
#Ӧ��avereps����������limma���������ݽ��д������������ͨ������������Ϣѧ�У�
#���ж��̽����߱�ǩ��Ӧͬһ������ʱ�������Ժϲ���Щ�в�����ƽ��ֵ��
data=avereps(data)
#������data��ƽ��ֵ����0���У�ɸѡ����һ���������Ļ���
data=data[rowMeans(data)>0,]

#��ȡĿ¼������"s1.txt"��β���ļ�
sampleName1=c()
files=dir()
files=grep("s1.txt$", files, value=T)
for(file in files){
    rt=read.table(file, header=F, sep="\t", check.names=F)      #��ȡ�����ļ�
    geneNames=as.vector(rt[,1])      #��ȡ��������
    uniqGene=unique(geneNames)       #����ȡunique
    sampleName1=c(sampleName1, uniqGene)
}

#��ȡĿ¼������"s2.txt"��β���ļ�
sampleName2=c()
files=dir()
files=grep("s2.txt$", files, value=T)
for(file in files){
    rt=read.table(file, header=F, sep="\t", check.names=F)      #��ȡ�����ļ�
    geneNames=as.vector(rt[,1])      #��ȡ��������
    uniqGene=unique(geneNames)       #����ȡunique
    sampleName2=c(sampleName2, uniqGene)
}

#��ȡʵ����Ͷ����������
conData=data[,sampleName1]
treatData=data[,sampleName2]
data=cbind(conData,treatData)
conNum=ncol(conData)
treatNum=ncol(treatData)

#�������
#������һ���ַ�����������"con"�ظ�conNum�Σ����������飻��"treat"�ظ�treatNum�Σ����������顣
Type=c(rep("con",conNum),rep("treat",treatNum))
#������һ����ƾ��󡣸þ�����������ģ����ϡ�~0��ʾģ����û�нؾ��factor(Type)���ַ�����Typeת��Ϊ���Ӳ���ģ��������Ԥ�������
design <- model.matrix(~0+factor(Type))
#Ϊ��ƾ����������
colnames(design) <- c("con","treat")
#ʹ��limma����lmFit������data��������������ݣ�����ƾ������������ϡ�
fit <- lmFit(data,design)
#�����ԱȾ������ڱȽ�"treat"��"con"�Ĳ���
cont.matrix<-makeContrasts(treat-con,levels=design)
#���ԱȾ���Ӧ������ǰ��ϵ�ģ�͡�
fit2 <- contrasts.fit(fit, cont.matrix)
#Ӧ�þ���Bayes��������ͳ�Ƶ���������limma�����еķ�������������߲�����������ͳ�ƹ�Ч��
fit2 <- eBayes(fit2)
#����Ͻ������ȡ����Ĳ�����ͳ����Ϣ�������˼ٷ����ʣ�FDR��У������ѡ�������200,000���������Ļ���
allDiff=topTable(fit2,adjust='fdr',number=200000)
#�������ӵ�allDiff���ݿ�Ķ������Դ���һ���µ����ݿ�allDiffOut��
allDiffOut=rbind(id=colnames(allDiff),allDiff)
#��allDiffOut���ݿ�д����Ϊ"all.txt"���ļ���ʹ���Ʊ�����Ϊ�ָ��������������ú�������
write.table(allDiffOut, file="all.txt", sep="\t", quote=F, col.names=F)

#���������ı�����
outData=rbind(id=paste0(colnames(data),"_",Type),data)
write.table(outData, file="normalize.txt", sep="\t", quote=F, col.names=F)

#���������
diffSig=allDiff[with(allDiff, (abs(logFC)>logFCfilter & adj.P.Val < adj.P.Val.Filter )), ]
diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut, file="diff.txt", sep="\t", quote=F, col.names=F)

#���������������
diffGeneExp=data[row.names(diffSig),]
diffGeneExpOut=rbind(id=paste0(colnames(diffGeneExp),"_",Type), diffGeneExp)
write.table(diffGeneExpOut, file="diffGeneExp.txt", sep="\t", quote=F, col.names=F)

#���Ʋ��������ͼ
geneNum=50
diffSig=diffSig[order(as.numeric(as.vector(diffSig$logFC))),]
#��������diffSig���ݿ򣬸���logFC�е�ֵ��
diffGeneName=as.vector(rownames(diffSig))
#��ȡ������diffSig���ݿ������������������������Ϊ����diffGeneName��
diffLength=length(diffGeneName)
#����diffGeneName�����ĳ��ȣ�Ҳ���������������������䱣����diffLength�С�
hmGene=c()
#��ʼ��һ���յ�����hmGene�������ڱ��潫Ҫ����ͼ����ʾ�Ļ�������
if(diffLength>(2*geneNum)){
    hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
}else{
    hmGene=diffGeneName
}
#�������������������������geneNum��Ҳ����100������ô�Ӳ��������ѡ��ǰ50���ͺ�50��������ѡ�����в������
hmExp=data[hmGene,]
#��ԭʼ��������ȡ��Ҫ����ͼ��չʾ�Ļ���ı���ֵ����������hmExp�С�
Type=c(rep("Con",conNum),rep("Treat",treatNum))
names(Type)=colnames(data)
#ΪType�����������֣���Щ��������data��������
Type=as.data.frame(Type)
#��Type����ת��Ϊ���ݿ�
pdf(file="heatmap.pdf", width=10, height=8)
#�������ͼ���豸ΪPDF����ָ������ļ����������ÿ��Ⱥ͸߶ȡ�
#ʹ��pheatmap����������ͼ��
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 8,
         fontsize_row=7,
         fontsize_col=8)
#hmExp��Ҫ��ʾ������
#annotation=Type��Ϊÿһ������ע�ͣ���"Con"��"Treat"��
#color��������ɫ����������ɫ����ɫ���ٵ���ɫ
#cluster_cols=F�������н��о���
#show_colnames=F������ʾ����
#scale="row"�����н������ţ�ʹÿһ�е�ֵ������ͬ�ľ�ֵ�ͱ�׼��
#fontsize��fontsize_row��fontsize_col�����������С		 
dev.off()
#�رյ�ǰ��ͼ���豸