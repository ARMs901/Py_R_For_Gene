# ����pROC�������û�а�װ����Ҫ�Ȱ�װ
# install.packages("pROC")
library(pROC) 

# ������������ļ��ͽ��������б��ļ���·��
expFile="diffGeneExp.txt" 
geneFile="intersected_genes.txt"

# ���ù���Ŀ¼
setwd("D:\\BaiduNetdiskDownload\\Diagnostic_video\\136Diagnostic\\16.ROC\\U_D")

# �ӱ��������ļ���ȡ����
rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)

# ����������ȡ�������ͣ����磺"con"����������ת��Ϊ������ֵ��0��1��
y=gsub("(.*)\\_(.*)", "\\2", colnames(rt))
y=ifelse(y=="con", 0, 1)

# �ӽ��������ļ��ж�ȡ�����б�
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)

# �������������б���Ϊÿ���������ROC����
for(x in as.vector(geneRT[,1])){
    # ʹ��y���������ͣ���rt[x,]������������ݣ�����ROC
    roc1=roc(y, as.numeric(rt[x,]))
    
    # ����ROC�����������AUC������������
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    
    # Ϊ��ǰ���򱣴�ROC���ߵ�PDF�ļ�
    pdf(file=paste0("ROC.",x,".pdf"), width=5, height=5)
    
    # ����ROC���ߣ�����ͼ����ʾAUCֵ
    plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=x)
    text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
    
    # �ر�PDF�豸
    dev.off()
}