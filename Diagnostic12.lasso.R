
#install.packages("glmnet")

set.seed(123) # �������������ȷ������Ŀ��ظ���

library(glmnet) # ����glmnet�����ð��ṩLasso�ع鹦��

inputFile="diffGeneExp.txt" # ָ�������ļ���

# ���ù���Ŀ¼����Ŀ¼Ӧ���������ļ�
setwd("D:\\BaiduNetdiskDownload\\Diagnostic_video\\136Diagnostic\\12.lasso\\U_D")

# ��ȡ�����ļ����ļ��ĵ�һ�б���Ϊ����������ת����ʹ����Ϊ�У�����Ϊ��
rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
rt=t(rt)

# ΪLassoģ�͹����������x����Ӧ����y
x=as.matrix(rt) # �����ݿ�ת��Ϊ����
y=gsub("(.*)\\_(.*)", "\\2", row.names(rt)) # ����������ȡ��Ӧ����

# ʹ��Lasso�����߼��ع�
fit=glmnet(x, y, family = "binomial", alpha=1)

# ʹ��10�۽�����֤ȷ����ѵ�lambdaֵ
cvfit=cv.glmnet(x, y, family="binomial", alpha=1,type.measure='deviance',nfolds = 10)

# ���������֤�Ľ��ΪPDFͼ��
pdf(file="cvfit.pdf",width=6,height=5.5)
plot(cvfit) # ����������֤�Ľ��
dev.off() # �ر�ͼ���豸

# ��ȡģ��ϵ������ȷ������ϵ����Ӧ������������Lassoѡ��Ļ���
coef=coef(fit, s = cvfit$lambda.min) # ��ȡ�����lambda�µ�ģ��ϵ��
index=which(coef != 0) # �ҳ�����ϵ��������
lassoGene=row.names(coef)[index] # ��ȡ��Ӧ�Ļ�����
lassoGene=lassoGene[-1] # ���б���ɾ���ؾ���

# ��ɸѡ���Ļ���д���ļ�
write.table(lassoGene, file="LASSO.gene.txt", sep="\t", quote=F, row.names=F, col.names=F)