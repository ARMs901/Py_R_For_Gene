# 检查是否已安装BiocManager包，若没有则进行安装
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# 使用BiocManager来安装limma包
#BiocManager::install("limma")
# 使用BiocManager来安装sva包
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("sva")
# 引用包
library(limma)  												# 引用limma包，用于生物统计分析
library(sva)    												# 引用sva包，用于批处理效应的校正
outFile="merge.txt"  											# 设定输出文件名
setwd("C:\\biowolf\\Diagnostic\\05.sva") 		 				# 设置工作目录

# 获取目录下所有".txt"结尾的文件
files = dir()  													# 获取当前目录下的所有文件
files = grep("txt$", files, value=T)  							# 筛选出.txt结尾的文件
geneList = list()  												# 初始化空列表用于存放基因信息

# 读取所有txt文件中的基因信息，保存到geneList
for(file in files){
    if(file == outFile){next}  									# 如果是输出文件，则跳过此次循环
    rt = read.table(file, header=T, sep="\t", check.names=F)  	# 读取文件内容
    geneNames = as.vector(rt[,1])  								# 获取基因名称
    uniqGene = unique(geneNames)   								# 对基因名称进行去重
    header = unlist(strsplit(file, "\\.|\\-"))  				# 分割文件名
    geneList[[header[1]]] = uniqGene  							# 将基因名称保存到列表
}
# 获取所有基因列表的交集
interGenes = Reduce(intersect, geneList)

# 初始化数据框和向量用于后续数据处理
allTab = data.frame()
batchType = c()
for(i in 1:length(files)){
    inputFile = files[i]
    header = unlist(strsplit(inputFile, "\\.|\\-"))  				# 分割文件名
    rt = read.table(inputFile, header=T, sep="\t", check.names=F)  	# 读取文件内容
    rt = as.matrix(rt)  											# 转换为矩阵格式
    rownames(rt) = rt[,1]  											# 设定矩阵的行名
    exp = rt[,2:ncol(rt)]  											# 获取除基因名称外的所有数据
    dimnames = list(rownames(exp), colnames(exp))  					# 获取行名和列名
    data = matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)  # 将数据转为数值矩阵
    rt = avereps(data)  											# 对数据进行处理

    # 对数值较大的数据进行log2变换
    qx = as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC = ( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
    if(LogC){
        rt[rt<0] = 0
        rt = log2(rt+1)
    }
    rt = normalizeBetweenArrays(rt)  								# 对数据进行归一化处理

    # 合并数据
    if(i == 1){
        allTab = rt[interGenes,]
    } else {
        allTab = cbind(allTab, rt[interGenes,])
    }
    batchType = c(batchType, rep(i, ncol(rt)))
}

# 对数据进行批处理效应的校正
outTab = ComBat(allTab, batchType, par.prior=TRUE)
outTab = rbind(geneNames=colnames(outTab), outTab)  					# 在数据的顶部添加列名
write.table(outTab, file="merge.txt", sep="\t", quote=F, col.names=F)  	# 将结果写入文件
