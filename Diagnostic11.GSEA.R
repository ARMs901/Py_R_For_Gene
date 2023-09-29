#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")


#���ð�
library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

inputFile="all.txt"         								   #�����ļ�
gmtFile="c2.cp.kegg.v7.4.symbols.gmt"      					   #�����ļ�
setwd("C:\\biowolf\\Diagnostic\\11.GSEA")      				   #���ù���Ŀ¼

# ��ȡ�����ļ��������������
rt = read.table(inputFile, header=T, sep="\t", check.names=F)  # ��ȡ�����ļ�Ϊ���ݿ�
logFC = as.vector(rt[,2])                                      # ��ȡ�ڶ��У�logFCֵ��Ϊ����
names(logFC) = as.vector(rt[,1])                               # ����һ������ΪlogFC����������
logFC = sort(logFC, decreasing=T)                              # ��logFCֵ���н�������

# ��������ļ�
gmt = read.gmt(gmtFile)                                        # ʹ��read.gmt������ȡgmt�ļ�

# ���и�������
kk = GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)              # ��logFC��������GSEA��������
kkTab = as.data.frame(kk)                                      # �����ת��Ϊ���ݿ�
kkTab = kkTab[kkTab$p.adjust<0.05,]                            # ѡ����Щ�������pֵС��0.05����
write.table(kkTab, file="GSEA.result.txt", sep="\t", quote=F, row.names = F)  # ��ɸѡ��Ľ��д��txt�ļ�

# ���ʵ���鸻����ͼ��
termNum = 5                                                    # ����Ҫչʾ��ͨ·��ĿΪ5
kkUp = kkTab[kkTab$NES>0,]                                     # ѡ�����Ĺ�һ�������÷֣�������ʵ�����и�����
if(nrow(kkUp) >= termNum){                                     # ��������ͨ·��Ŀ���ڻ����5
    showTerm = row.names(kkUp)[1:termNum]                      # ѡ��ǰ5��ͨ·
    gseaplot = gseaplot2(kk, showTerm, base_size=8, title="Enriched in Treat")  # ʹ��gseaplot2���Ƹ���ͼ
    pdf(file="GSEA.treat.pdf", width=7, height=5.5)            # �������ΪPDF��ʽ
    print(gseaplot)                                            # ���ͼ��
    dev.off()                                                  # �ر�ͼ���豸
}

# ��������鸻����ͼ��
termNum = 5                                                    # ����Ҫչʾ��ͨ·��ĿΪ5
kkDown = kkTab[kkTab$NES<0,]                                   # ѡ�񸺵Ĺ�һ�������÷֣��������������и�����
if(nrow(kkDown) >= termNum){                                   # ��������ͨ·��Ŀ���ڻ����5
    showTerm = row.names(kkDown)[1:termNum]                    # ѡ��ǰ5��ͨ·
    gseaplot = gseaplot2(kk, showTerm, base_size=8, title="Enriched in Control")  # ʹ��gseaplot2���Ƹ���ͼ
    pdf(file="GSEA.con.pdf", width=7, height=5.5)              # �������ΪPDF��ʽ
    print(gseaplot)                                            # ���ͼ��
    dev.off()                                                  # �ر�ͼ���豸
}