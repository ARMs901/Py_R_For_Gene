
#install.packages("e1071")
#install.packages("kernlab")
#install.packages("caret")

# 加载所需的R包
library(e1071)      # 加载e1071包，包含支持向量机和相关函数
library(kernlab)    # 加载kernlab包，用于核方法（包括SVM）
library(caret)      # 加载caret包，用于机器学习和数据建模

set.seed(123)       # 设置随机种子，确保结果的可重复性

inputFile="diffGeneExp.txt"  # 指定输入文件名

# 设置工作目录，该目录应包含输入文件
setwd("D:\\BaiduNetdiskDownload\\Diagnostic_video\\136Diagnostic\\13.SVM\\U_D")

# 读取输入文件。文件的第一列被视为行名，数据转置以使基因为列，样本为行
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(data)

# 从行名中提取样本的类别标签
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))

# 使用SVM-RFE进行特征选择
Profile=rfe(x=data,
            y=as.numeric(as.factor(group)),  # 将类别标签转化为数值格式
            sizes = c(2,4,6,8, seq(10,40,by=3)),  # 设定考虑的特征数量
            rfeControl = rfeControl(functions = caretFuncs, method = "cv"),  # 使用交叉验证评估模型性能
            methods="svmRadial")  # 使用径向基函数核的SVM

# 输出SVM-RFE结果的图形
pdf(file="SVM-RFE.pdf", width=6, height=5.5)  # 设置输出为PDF格式
par(las=1)
x = Profile$results$Variables  # 提取考虑的特征数量
y = Profile$results$RMSE       # 提取交叉验证的均方根误差
plot(x, y, xlab="Variables", ylab="RMSE (Cross-Validation)", col="darkgreen")  # 绘制RMSE随特征数量的变化曲线
lines(x, y, col="darkgreen")
# 标注交叉验证误差最小的点
wmin=which.min(y)  # 找到RMSE最小的索引
wmin.x=x[wmin]
wmin.y=y[wmin]
points(wmin.x, wmin.y, col="blue", pch=16)  # 标注点
text(wmin.x, wmin.y, paste0('N=',wmin.x), pos=2, col=2)  # 显示考虑的特征数量
dev.off()  # 关闭图形设备

# 输出被选择的特征（基因）
featureGenes=Profile$optVariables  # 提取最佳的特征集
write.table(file="SVM-RFE.gene.txt", featureGenes, sep="\t", quote=F, row.names=F, col.names=F)  # 将特征写入txt文件
