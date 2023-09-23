# ����Ƿ��Ѱ�װBiocManager������û������а�װ
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# ʹ��BiocManager����װlimma��
#BiocManager::install("limma")
# ʹ��BiocManager����װsva��
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("sva")
# ���ð�
library(limma)  												# ����limma������������ͳ�Ʒ���
library(sva)    												# ����sva��������������ЧӦ��У��
outFile="merge.txt"  											# �趨����ļ���
setwd("C:\\biowolf\\Diagnostic\\05.sva") 		 				# ���ù���Ŀ¼

# ��ȡĿ¼������".txt"��β���ļ�
files = dir()  													# ��ȡ��ǰĿ¼�µ������ļ�
files = grep("txt$", files, value=T)  							# ɸѡ��.txt��β���ļ�
geneList = list()  												# ��ʼ�����б����ڴ�Ż�����Ϣ

# ��ȡ����txt�ļ��еĻ�����Ϣ�����浽geneList
for(file in files){
    if(file == outFile){next}  									# ���������ļ����������˴�ѭ��
    rt = read.table(file, header=T, sep="\t", check.names=F)  	# ��ȡ�ļ�����
    geneNames = as.vector(rt[,1])  								# ��ȡ��������
    uniqGene = unique(geneNames)   								# �Ի������ƽ���ȥ��
    header = unlist(strsplit(file, "\\.|\\-"))  				# �ָ��ļ���
    geneList[[header[1]]] = uniqGene  							# ���������Ʊ��浽�б�
}
# ��ȡ���л����б��Ľ���
interGenes = Reduce(intersect, geneList)

# ��ʼ�����ݿ���������ں������ݴ���
allTab = data.frame()
batchType = c()
for(i in 1:length(files)){
    inputFile = files[i]
    header = unlist(strsplit(inputFile, "\\.|\\-"))  				# �ָ��ļ���
    rt = read.table(inputFile, header=T, sep="\t", check.names=F)  	# ��ȡ�ļ�����
    rt = as.matrix(rt)  											# ת��Ϊ�����ʽ
    rownames(rt) = rt[,1]  											# �趨���������
    exp = rt[,2:ncol(rt)]  											# ��ȡ���������������������
    dimnames = list(rownames(exp), colnames(exp))  					# ��ȡ����������
    data = matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)  # ������תΪ��ֵ����
    rt = avereps(data)  											# �����ݽ��д���

    # ����ֵ�ϴ�����ݽ���log2�任
    qx = as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC = ( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
    if(LogC){
        rt[rt<0] = 0
        rt = log2(rt+1)
    }
    rt = normalizeBetweenArrays(rt)  								# �����ݽ��й�һ������

    # �ϲ�����
    if(i == 1){
        allTab = rt[interGenes,]
    } else {
        allTab = cbind(allTab, rt[interGenes,])
    }
    batchType = c(batchType, rep(i, ncol(rt)))
}

# �����ݽ���������ЧӦ��У��
outTab = ComBat(allTab, batchType, par.prior=TRUE)
outTab = rbind(geneNames=colnames(outTab), outTab)  					# �����ݵĶ�����������
write.table(outTab, file="merge.txt", sep="\t", quote=F, col.names=F)  	# �����д���ļ�