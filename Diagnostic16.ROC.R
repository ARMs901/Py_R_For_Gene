# 引用pROC包，如果没有安装，需要先安装
# install.packages("pROC")
library(pROC) 

# 定义表达数据文件和交集基因列表文件的路径
expFile="diffGeneExp.txt" 
geneFile="intersected_genes.txt"

# 设置工作目录
setwd("D:\\BaiduNetdiskDownload\\Diagnostic_video\\136Diagnostic\\16.ROC\\U_D")

# 从表达数据文件读取数据
rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)

# 从列名中提取样本类型（例如："con"），并将其转换为二进制值（0或1）
y=gsub("(.*)\\_(.*)", "\\2", colnames(rt))
y=ifelse(y=="con", 0, 1)

# 从交集基因文件中读取基因列表
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)

# 遍历交集基因列表，为每个基因绘制ROC曲线
for(x in as.vector(geneRT[,1])){
    # 使用y（样本类型）和rt[x,]（基因表达数据）计算ROC
    roc1=roc(y, as.numeric(rt[x,]))
    
    # 计算ROC曲线下面积（AUC）的置信区间
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    
    # 为当前基因保存ROC曲线到PDF文件
    pdf(file=paste0("ROC.",x,".pdf"), width=5, height=5)
    
    # 绘制ROC曲线，并在图上显示AUC值
    plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=x)
    text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
    
    # 关闭PDF设备
    dev.off()
}
