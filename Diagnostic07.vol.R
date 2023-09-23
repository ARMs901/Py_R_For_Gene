
#install.packages("ggplot2")
#install.packages("ggrepel")

#���ð�
library(dplyr)
library(ggplot2)
library(ggrepel)

logFCfilter=2               #logFC��������
adj.P.Val.Filter=0.05       #�������pֵ��������
inputFile="all.txt"         #�����ļ�
setwd("C:\\biowolf\\Diagnostic\\07.vol")       #���ù���Ŀ¼

#��ȡ�����ļ�
rt = read.table(inputFile, header=T, sep="\t", check.names=F)
#���������ԣ�����adj.P.Val��logFC����ֵ����һ���±���Sig�����ڱ�ʾ�����Ƿ��ϵ����µ����ǲ�������
Sig=ifelse((rt$adj.P.Val<adj.P.Val.Filter) & (abs(rt$logFC)>logFCfilter), ifelse(rt$logFC>logFCfilter,"Up","Down"), "Not")

#���ƻ�ɽͼ
rt = mutate(rt, Sig=Sig)#ʹ��mutate����������dplyr�������´�����Sig�����ӵ����ݿ�rt�С�
#ʹ��ggplot2��������ɽͼ�����У�
#aes(logFC, -log10(adj.P.Val))���趨x��ΪlogFC��y��Ϊ-log10(adj.P.Val)��
#geom_point(aes(col=Sig))�����Ƶ㣬������Sig��Ϊ����ɫ��
#scale_color_manual���ֶ�Ϊ���趨��ɫ��
#labs��theme������ͼ�ı������ʽ��
p = ggplot(rt, aes(logFC, -log10(adj.P.Val)))+
    geom_point(aes(col=Sig))+
    scale_color_manual(values=c("green", "black","red"))+
    labs(title = " ")+
    theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
#�ڻ�ɽͼ��Ϊ�����Ļ������ӱ�ǩ��ʹ��geom_label_repel������ggrepel����ȷ����ǩ֮�䲻�ص���
#���ڲ��������Ļ��򣬱�ע���������
p1=p+geom_label_repel(data=filter(rt, ((rt$adj.P.Val<adj.P.Val.Filter) & (abs(rt$logFC)>logFCfilter))),
                    box.padding=0.1, point.padding=0.1, min.segment.length=0.05,
                    size=1.8, aes(label=id)) + theme_bw()
#�����ɽͼ
pdf(file="vol.pdf", width=7, height=6.1)
print(p1)
dev.off()

