import pandas as pd
import os

# 定义输入文件的路径
ann_file_path = "D:/BaiduNetdiskDownload/Diagnostic_video/136Diagnostic/04.ann/GPL19833/ann.xlsx"
probe_file_path = "D:/BaiduNetdiskDownload/Diagnostic_video/136Diagnostic/04.ann/GPL19833/probe.xlsx"

# 读取输入文件
ann_df = pd.read_excel(ann_file_path, usecols=["ID", "GENE_SYMBOL"])
probe_df = pd.read_excel(probe_file_path)

# 合并两个数据框
df = pd.merge(ann_df, probe_df, how="inner", left_on="ID", right_on="ID_REF")

# 删除没有基因注释的探针
df = df[df["GENE_SYMBOL"].notna()]

# 删除 ID 和 ID_REF 这两列
df = df.drop(["ID", "ID_REF"], axis=1)

# 输出结果
df.to_csv(os.path.join(os.path.dirname(ann_file_path), "result.csv"), index=False)
