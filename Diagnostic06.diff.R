
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
#install.packages("pheatmap")

#引用包
library(limma)
library(pheatmap)

inputFile="merge.txt"       #输入文件
logFCfilter=2               #logFC过滤阈值
adj.P.Val.Filter=0.05       #矫正后p值阈值
setwd("C:\\biowolf\\Diagnostic\\06.diff")      #设置工作目录

#读取输入文件，并对输入文件整理
#使用read.table函数从inputFile中读取数据。其中：
#header=T：文件的第一行是列名。
#sep="\t"：字段是用制表符（tab）分隔的。
#check.names=F：不修改列名以使其变为标准的R变量名（例如，如果列名包含空格或特殊字符，R通常会修改它们）。
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
#将rt转换为一个矩阵。
rt=as.matrix(rt)
#将rt矩阵的第一列设置为行名。
rownames(rt)=rt[,1]
#提取rt矩阵的第二列到最后一列，并将结果保存在exp中。
exp=rt[,2:ncol(rt)]
#为下一步的数据矩阵创建维度名称。这里获取exp的行名和列名，并将它们放入一个列表中，保存在dimnames中。
dimnames=list(rownames(exp),colnames(exp))
#将exp再次转换为一个数值矩阵（确保所有元素都是数字），并为这个新矩阵设置行名和列名。
#这样做是为了确保所有的元素都是数值类型。
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
#应用avereps函数（来自limma包）对数据进行处理。这个函数通常用于生物信息学中，
#当有多个探针或者标签对应同一个基因时，它可以合并这些行并计算平均值。
data=avereps(data)
#仅保留data中平均值大于0的行，筛选出有一定表达量的基因。
data=data[rowMeans(data)>0,]

#读取目录下所有"s1.txt"结尾的文件
sampleName1=c()
files=dir()
files=grep("s1.txt$", files, value=T)
for(file in files){
    rt=read.table(file, header=F, sep="\t", check.names=F)      #读取输入文件
    geneNames=as.vector(rt[,1])      #提取基因名称
    uniqGene=unique(geneNames)       #基因取unique
    sampleName1=c(sampleName1, uniqGene)
}

#读取目录下所有"s2.txt"结尾的文件
sampleName2=c()
files=dir()
files=grep("s2.txt$", files, value=T)
for(file in files){
    rt=read.table(file, header=F, sep="\t", check.names=F)      #读取输入文件
    geneNames=as.vector(rt[,1])      #提取基因名称
    uniqGene=unique(geneNames)       #基因取unique
    sampleName2=c(sampleName2, uniqGene)
}

#提取实验组和对照组的数据
valid_sampleName1 <- sampleName1[sampleName1 %in% colnames(data)]
valid_sampleName2 <- sampleName2[sampleName2 %in% colnames(data)]
#让不存在与merge.txt里的samplename进行舍弃，让存在的继续运行下去，不让代码报错
# 现在，使用这些有效的列名提取数据
conData <- data[, valid_sampleName1]
treatData <- data[, valid_sampleName2]
data=cbind(conData,treatData)
conNum=ncol(conData)
treatNum=ncol(treatData)

#差异分析
#创建了一个字符向量，其中"con"重复conNum次，代表对照组；而"treat"重复treatNum次，代表处理组。
Type=c(rep("con",conNum),rep("treat",treatNum))
#建立了一个设计矩阵。该矩阵将用于线性模型拟合。~0表示模型中没有截距项。factor(Type)将字符向量Type转化为因子并在模型中用作预测变量。
design <- model.matrix(~0+factor(Type))
#为设计矩阵的列命名
colnames(design) <- c("con","treat")
#使用limma包的lmFit函数将data（即基因表达数据）与设计矩阵进行线性拟合。
fit <- lmFit(data,design)
#建立对比矩阵，用于比较"treat"和"con"的差异
cont.matrix<-makeContrasts(treat-con,levels=design)
#将对比矩阵应用于先前拟合的模型。
fit2 <- contrasts.fit(fit, cont.matrix)
#应用经验Bayes方法进行统计调整。这是limma包特有的方法，它可以提高差异表达分析的统计功效。
fit2 <- eBayes(fit2)
#从拟合结果中提取表达的差异性统计信息。考虑了假发现率（FDR）校正，并选择了最多200,000个差异表达的基因。
allDiff=topTable(fit2,adjust='fdr',number=200000)
#将列名加到allDiff数据框的顶部，以创建一个新的数据框allDiffOut。
allDiffOut=rbind(id=colnames(allDiff),allDiff)
#将allDiffOut数据框写入名为"all.txt"的文件，使用制表符作为分隔符，并避免引用和列名。
write.table(allDiffOut, file="all.txt", sep="\t", quote=F, col.names=F)

#输出矫正后的表达量
outData=rbind(id=paste0(colnames(data),"_",Type),data)
write.table(outData, file="normalize.txt", sep="\t", quote=F, col.names=F)

#输出差异结果
diffSig=allDiff[with(allDiff, (abs(logFC)>logFCfilter & adj.P.Val < adj.P.Val.Filter )), ]
diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut, file="diff.txt", sep="\t", quote=F, col.names=F)

#输出差异基因表达量
diffGeneExp=data[row.names(diffSig),]
diffGeneExpOut=rbind(id=paste0(colnames(diffGeneExp),"_",Type), diffGeneExp)
write.table(diffGeneExpOut, file="diffGeneExp.txt", sep="\t", quote=F, col.names=F)

#绘制差异基因热图
geneNum=50
diffSig=diffSig[order(as.numeric(as.vector(diffSig$logFC))),]
#重新排序diffSig数据框，根据logFC列的值。
diffGeneName=as.vector(rownames(diffSig))
#提取排序后的diffSig数据框的行名（即基因名）并保存为向量diffGeneName。
diffLength=length(diffGeneName)
#计算diffGeneName向量的长度，也即差异基因的数量，并将其保存在diffLength中。
hmGene=c()
#初始化一个空的向量hmGene，将用于保存将要在热图中显示的基因名。
if(diffLength>(2*geneNum)){
    hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
}else{
    hmGene=diffGeneName
}
#如果差异基因的数量大于两倍的geneNum（也就是100），那么从差异基因中选择前50个和后50个。否则，选择所有差异基因。
hmExp=data[hmGene,]
#从原始数据中提取将要在热图中展示的基因的表达值，并保存在hmExp中。
Type=c(rep("Con",conNum),rep("Treat",treatNum))
names(Type)=colnames(data)
#为Type向量设置名字，这些名字来自data的列名。
Type=as.data.frame(Type)
#将Type向量转换为数据框。
pdf(file="heatmap.pdf", width=10, height=8)
#设置输出图形设备为PDF，并指定输出文件名，并设置宽度和高度。
#使用pheatmap函数生成热图。
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 8,
         fontsize_row=7,
         fontsize_col=8)
#hmExp：要显示的数据
#annotation=Type：为每一列添加注释（即"Con"或"Treat"）
#color：设置颜色方案，从蓝色到白色，再到红色
#cluster_cols=F：不对列进行聚类
#show_colnames=F：不显示列名
#scale="row"：按行进行缩放，使每一行的值具有相同的均值和标准差
#fontsize、fontsize_row和fontsize_col：设置字体大小		 
dev.off()
#关闭当前的图形设备
