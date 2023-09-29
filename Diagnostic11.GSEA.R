#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")

#引用包
library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

inputFile="all.txt"         								   #输入文件
gmtFile="c2.cp.kegg.v7.4.symbols.gmt"      					   #基因集文件
setwd("C:\\biowolf\\Diagnostic\\11.GSEA")      				   #设置工作目录

# 读取输入文件并对其进行整理
rt = read.table(inputFile, header=T, sep="\t", check.names=F)  # 读取输入文件为数据框
logFC = as.vector(rt[,2])                                      # 提取第二列（logFC值）为向量
names(logFC) = as.vector(rt[,1])                               # 将第一列设置为logFC向量的名字
logFC = sort(logFC, decreasing=T)                              # 对logFC值进行降序排序

# 读入基因集文件
gmt = read.gmt(gmtFile)                                        # 使用read.gmt函数读取gmt文件

# 进行富集分析
kk = GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)              # 对logFC向量进行GSEA富集分析
kkTab = as.data.frame(kk)                                      # 将结果转化为数据框
kkTab = kkTab[kkTab$p.adjust<0.05,]                            # 选择那些调整后的p值小于0.05的行
write.table(kkTab, file="GSEA.result.txt", sep="\t", quote=F, row.names = F)  # 将筛选后的结果写入txt文件

# 输出实验组富集的图形
termNum = 5                                                    # 设置要展示的通路数目为5
kkUp = kkTab[kkTab$NES>0,]                                     # 选择正的归一化富集得分（代表在实验组中富集）
if(nrow(kkUp) >= termNum){                                     # 如果满足的通路数目大于或等于5
    showTerm = row.names(kkUp)[1:termNum]                      # 选择前5个通路
    gseaplot = gseaplot2(kk, showTerm, base_size=8, title="Enriched in Treat")  # 使用gseaplot2绘制富集图
    pdf(file="GSEA.treat.pdf", width=7, height=5.5)            # 设置输出为PDF格式
    print(gseaplot)                                            # 输出图形
    dev.off()                                                  # 关闭图形设备
}

# 输出正常组富集的图形
termNum = 5                                                    # 设置要展示的通路数目为5
kkDown = kkTab[kkTab$NES<0,]                                   # 选择负的归一化富集得分（代表在正常组中富集）
if(nrow(kkDown) >= termNum){                                   # 如果满足的通路数目大于或等于5
    showTerm = row.names(kkDown)[1:termNum]                    # 选择前5个通路
    gseaplot = gseaplot2(kk, showTerm, base_size=8, title="Enriched in Control")  # 使用gseaplot2绘制富集图
    pdf(file="GSEA.con.pdf", width=7, height=5.5)              # 设置输出为PDF格式
    print(gseaplot)                                            # 输出图形
    dev.off()                                                  # 关闭图形设备
}
