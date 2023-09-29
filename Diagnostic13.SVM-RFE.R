
#install.packages("e1071")
#install.packages("kernlab")
#install.packages("caret")

# ���������R��
library(e1071)      # ����e1071��������֧������������غ���
library(kernlab)    # ����kernlab�������ں˷���������SVM��
library(caret)      # ����caret�������ڻ���ѧϰ�����ݽ�ģ

set.seed(123)       # ����������ӣ�ȷ������Ŀ��ظ���

inputFile="diffGeneExp.txt"  # ָ�������ļ���

# ���ù���Ŀ¼����Ŀ¼Ӧ���������ļ�
setwd("D:\\BaiduNetdiskDownload\\Diagnostic_video\\136Diagnostic\\13.SVM\\U_D")

# ��ȡ�����ļ����ļ��ĵ�һ�б���Ϊ����������ת����ʹ����Ϊ�У�����Ϊ��
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(data)

# ����������ȡ����������ǩ
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))

# ʹ��SVM-RFE��������ѡ��
Profile=rfe(x=data,
            y=as.numeric(as.factor(group)),  # ������ǩת��Ϊ��ֵ��ʽ
            sizes = c(2,4,6,8, seq(10,40,by=3)),  # �趨���ǵ���������
            rfeControl = rfeControl(functions = caretFuncs, method = "cv"),  # ʹ�ý�����֤����ģ������
            methods="svmRadial")  # ʹ�þ���������˵�SVM

# ���SVM-RFE�����ͼ��
pdf(file="SVM-RFE.pdf", width=6, height=5.5)  # �������ΪPDF��ʽ
par(las=1)
x = Profile$results$Variables  # ��ȡ���ǵ���������
y = Profile$results$RMSE       # ��ȡ������֤�ľ��������
plot(x, y, xlab="Variables", ylab="RMSE (Cross-Validation)", col="darkgreen")  # ����RMSE�����������ı仯����
lines(x, y, col="darkgreen")
# ��ע������֤�����С�ĵ�
wmin=which.min(y)  # �ҵ�RMSE��С������
wmin.x=x[wmin]
wmin.y=y[wmin]
points(wmin.x, wmin.y, col="blue", pch=16)  # ��ע��
text(wmin.x, wmin.y, paste0('N=',wmin.x), pos=2, col=2)  # ��ʾ���ǵ���������
dev.off()  # �ر�ͼ���豸

# �����ѡ�������������
featureGenes=Profile$optVariables  # ��ȡ��ѵ�������
write.table(file="SVM-RFE.gene.txt", featureGenes, sep="\t", quote=F, row.names=F, col.names=F)  # ������д��txt�ļ�