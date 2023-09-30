# 加载 vioplot 包，用于绘制小提琴图
library(vioplot)

# 定义输入文件的路径
inputFile="CIBERSORT-Results.txt"

# 设置工作目录
setwd("D:\\BaiduNetdiskDownload\\Diagnostic_video\\136Diagnostic\\20.vioplot\\U_D")

# 从文件中读取免疫细胞浸润数据
rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)

# 使用正则表达式对样品进行分组：控制组 (con) 和治疗组 (treat)
con=grepl("_con", rownames(rt), ignore.case=T)
treat=grepl("_treat", rownames(rt), ignore.case=T)
conData=rt[con,]
treatData=rt[treat,]
conNum=nrow(conData)
treatNum=nrow(treatData)
rt=rbind(conData,treatData)

# 初始化输出数据框
outTab=data.frame()

# 设置输出 PDF 文件的尺寸
pdf(file="vioplot.pdf", height=8, width=13)
par(las=1,mar=c(10,6,3,3))

# 创建一个空白的图形
x=c(1:ncol(rt))
y=c(1:ncol(rt))
plot(x,y,
     xlim=c(0,63),ylim=c(min(rt),max(rt)+0.05),
     main="",xlab="", ylab="Fraction",
     pch=21,
     col="white",
     xaxt="n")

# 对每个免疫细胞循环，绘制小提琴图，对照组用蓝色表示，实验组用红色表示
for(i in 1:ncol(rt)){
  # 检查标准差是否为零，避免绘图错误
  if(sd(rt[1:conNum,i])==0){
    rt[1,i]=0.00001
  }
  if(sd(rt[(conNum+1):(conNum+treatNum),i])==0){
    rt[(conNum+1),i]=0.00001
  }

  # 分别提取控制组和治疗组的数据
  conData=rt[1:conNum,i]
  treatData=rt[(conNum+1):(conNum+treatNum),i]

  # 绘制小提琴图
  vioplot(conData,at=3*(i-1),lty=1,add = T,col = 'blue')
  vioplot(treatData,at=3*(i-1)+1,lty=1,add = T,col = 'red')

  # 进行威尔科克森秩和检验
  wilcoxTest=wilcox.test(conData,treatData)
  p=wilcoxTest$p.value
  
  # 如果 p 值小于 0.05，则保存到输出数据框
  if(p<0.05){
    cellPvalue=cbind(Cell=colnames(rt)[i],pvalue=p)
    outTab=rbind(outTab,cellPvalue)
  }

  # 在小提琴图上添加 p 值标签
  mx=max(c(conData,treatData))
  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
  text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",sprintf("%.03f",p))), cex = 0.8)
}

# 添加图例
legend("topright", 
       c("Con", "Treat"),
       lwd=3,bty="n",cex=1,
       col=c("blue","red"))

# 添加 x 轴标签
text(seq(1,64,3),-0.05,xpd = NA,labels=colnames(rt),cex = 1,srt = 45,pos=2)

# 关闭图形输出
dev.off()

# 将显著差异的免疫细胞和它们的 p 值保存到文件
write.table(outTab,file="immuneDiff.xls",sep="\t",row.names=F,quote=F)
