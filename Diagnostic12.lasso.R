
#install.packages("glmnet")

set.seed(123) # 设置随机种子以确保结果的可重复性

library(glmnet) # 加载glmnet包，该包提供Lasso回归功能

inputFile="diffGeneExp.txt" # 指定输入文件名

# 设置工作目录，该目录应包含输入文件
setwd("D:\\BaiduNetdiskDownload\\Diagnostic_video\\136Diagnostic\\12.lasso\\U_D")

# 读取输入文件。文件的第一列被视为行名，数据转置以使基因为列，样本为行
rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
rt=t(rt)

# 为Lasso模型构建输入矩阵x和响应变量y
x=as.matrix(rt) # 将数据框转换为矩阵
y=gsub("(.*)\\_(.*)", "\\2", row.names(rt)) # 从行名中提取响应变量

# 使用Lasso进行逻辑回归
fit=glmnet(x, y, family = "binomial", alpha=1)

# 使用10折交叉验证确定最佳的lambda值
cvfit=cv.glmnet(x, y, family="binomial", alpha=1,type.measure='deviance',nfolds = 10)

# 输出交叉验证的结果为PDF图形
pdf(file="cvfit.pdf",width=6,height=5.5)
plot(cvfit) # 画出交叉验证的结果
dev.off() # 关闭图形设备

# 获取模型系数，并确定非零系数对应的特征（即被Lasso选择的基因）
coef=coef(fit, s = cvfit$lambda.min) # 获取在最佳lambda下的模型系数
index=which(coef != 0) # 找出非零系数的索引
lassoGene=row.names(coef)[index] # 获取对应的基因名
lassoGene=lassoGene[-1] # 从列表中删除截距项

# 将筛选出的基因写入文件
write.table(lassoGene, file="LASSO.gene.txt", sep="\t", quote=F, row.names=F, col.names=F)
