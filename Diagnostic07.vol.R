
#install.packages("ggplot2")
#install.packages("ggrepel")

#引用包
library(dplyr)
library(ggplot2)
library(ggrepel)

logFCfilter=2               #logFC过滤条件
adj.P.Val.Filter=0.05       #矫正后的p值过滤条件
inputFile="all.txt"         #输入文件
setwd("C:\\biowolf\\Diagnostic\\07.vol")       #设置工作目录

#读取输入文件
rt = read.table(inputFile, header=T, sep="\t", check.names=F)
#定义显著性：基于adj.P.Val和logFC的阈值创建一个新变量Sig，用于表示基因是否上调、下调还是不显著。
Sig=ifelse((rt$adj.P.Val<adj.P.Val.Filter) & (abs(rt$logFC)>logFCfilter), ifelse(rt$logFC>logFCfilter,"Up","Down"), "Not")

#绘制火山图
rt = mutate(rt, Sig=Sig)#使用mutate函数（来自dplyr包）将新创建的Sig列添加到数据框rt中。
#使用ggplot2包创建火山图。其中：
#aes(logFC, -log10(adj.P.Val))：设定x轴为logFC，y轴为-log10(adj.P.Val)。
#geom_point(aes(col=Sig))：绘制点，并根据Sig列为点上色。
#scale_color_manual：手动为点设定颜色。
#labs和theme：设置图的标题和样式。
p = ggplot(rt, aes(logFC, -log10(adj.P.Val)))+
    geom_point(aes(col=Sig))+
    scale_color_manual(values=c("green", "black","red"))+
    labs(title = " ")+
    theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
#在火山图上为显著的基因添加标签。使用geom_label_repel（来自ggrepel包）确保标签之间不重叠。
#对于差异显著的基因，标注基因的名称
p1=p+geom_label_repel(data=filter(rt, ((rt$adj.P.Val<adj.P.Val.Filter) & (abs(rt$logFC)>logFCfilter))),
                    box.padding=0.1, point.padding=0.1, min.segment.length=0.05,
                    size=1.8, aes(label=id)) + theme_bw()
#输出火山图
pdf(file="vol.pdf", width=7, height=6.1)
print(p1)
dev.off()


