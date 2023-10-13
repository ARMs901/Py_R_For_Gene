library(biomaRt)
library(readxl)
library(writexl)

# 选择mart
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# 读取Excel文件
file_path <- "D:/BaiduNetdiskDownload/Diagnostic_video/136Diagnostic/04.ann/UC/UC_GSE65114/ann.xlsx"
df <- read_xlsx(file_path)

# 使用biomaRt查询基因名称
results <- getBM(attributes = c('refseq_mrna', 'external_gene_name'), 
                 filters = 'refseq_mrna', 
                 values = df$GB_ACC, 
                 mart = mart)

# 将查询结果合并回原始数据框
df <- merge(df, results, by.x = "GB_ACC", by.y = "refseq_mrna", all.x = TRUE)

# 保存结果到新的Excel文件
output_path <- sub(".xlsx", "_converted.xlsx", file_path)
write_xlsx(df, output_path)
