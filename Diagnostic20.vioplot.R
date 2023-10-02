# ���� vioplot �������ڻ���С����ͼ
library(vioplot)

# ���������ļ���·��
inputFile="CIBERSORT-Results.txt"

# ���ù���Ŀ¼
setwd("D:\\BaiduNetdiskDownload\\Diagnostic_video\\136Diagnostic\\20.vioplot\\U_D")

# ���ļ��ж�ȡ����ϸ����������
rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)

# ʹ���������ʽ����Ʒ���з��飺������ (con) �������� (treat)
con=grepl("_con", rownames(rt), ignore.case=T)
treat=grepl("_treat", rownames(rt), ignore.case=T)
conData=rt[con,]
treatData=rt[treat,]
conNum=nrow(conData)
treatNum=nrow(treatData)
rt=rbind(conData,treatData)

# ��ʼ��������ݿ�
outTab=data.frame()

# ������� PDF �ļ��ĳߴ�
pdf(file="vioplot.pdf", height=8, width=13)
par(las=1,mar=c(10,6,3,3))

# ����һ���հ׵�ͼ��
x=c(1:ncol(rt))
y=c(1:ncol(rt))
plot(x,y,
     xlim=c(0,63),ylim=c(min(rt),max(rt)+0.05),
     main="",xlab="", ylab="Fraction",
     pch=21,
     col="white",
     xaxt="n")

# ��ÿ������ϸ��ѭ��������С����ͼ������������ɫ��ʾ��ʵ�����ú�ɫ��ʾ
for(i in 1:ncol(rt)){
  # ����׼���Ƿ�Ϊ�㣬�����ͼ����
  if(sd(rt[1:conNum,i])==0){
    rt[1,i]=0.00001
  }
  if(sd(rt[(conNum+1):(conNum+treatNum),i])==0){
    rt[(conNum+1),i]=0.00001
  }

  # �ֱ���ȡ������������������
  conData=rt[1:conNum,i]
  treatData=rt[(conNum+1):(conNum+treatNum),i]

  # ����С����ͼ
  vioplot(conData,at=3*(i-1),lty=1,add = T,col = 'blue')
  vioplot(treatData,at=3*(i-1)+1,lty=1,add = T,col = 'red')

  # ���������ƿ�ɭ�Ⱥͼ���
  wilcoxTest=wilcox.test(conData,treatData)
  p=wilcoxTest$p.value
  
  # ��� p ֵС�� 0.05���򱣴浽������ݿ�
  if(p<0.05){
    cellPvalue=cbind(Cell=colnames(rt)[i],pvalue=p)
    outTab=rbind(outTab,cellPvalue)
  }

  # ��С����ͼ������ p ֵ��ǩ
  mx=max(c(conData,treatData))
  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
  text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",sprintf("%.03f",p))), cex = 0.8)
}

# ����ͼ��
legend("topright", 
       c("Con", "Treat"),
       lwd=3,bty="n",cex=1,
       col=c("blue","red"))

# ���� x ���ǩ
text(seq(1,64,3),-0.05,xpd = NA,labels=colnames(rt),cex = 1,srt = 45,pos=2)

# �ر�ͼ�����
dev.off()

# ���������������ϸ�������ǵ� p ֵ���浽�ļ�
write.table(outTab,file="immuneDiff.xls",sep="\t",row.names=F,quote=F)