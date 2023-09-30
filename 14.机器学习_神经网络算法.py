# 导入所需的库
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.neural_network import MLPClassifier

# 加载数据
# 使用 pandas 的 read_csv 方法从指定的路径加载数据文件。这个文件似乎是一个以 tab 分隔的文本文件，其中第一列作为索引。
data = pd.read_csv("D:/BaiduNetdiskDownload/Diagnostic_video/136Diagnostic/13.SVM/U_D/py_nn_feature/diffGeneExp.txt", sep="\t", index_col=0)
# 将数据转置，这样样本将作为行而基因作为列。
data_transposed = data.transpose()

# 从数据的行名 (样本名) 中提取样本类型，如 "con"。
labels = data_transposed.index.str.split("_").str[-1]
# 将样本标签转换为数字形式，其中 "con" 标签为 0，其他标签为 1。
numeric_labels = np.where(labels == "con", 0, 1)

# 初始化一个多层感知器 (神经网络) 分类器。
# 定义两个隐藏层，第一个有 50 个节点，第二个有 25 个节点。
# max_iter 设为 1000，表示最大迭代次数。
# random_state 设为 42 以确保每次运行都获得相同的结果。
mlp = MLPClassifier(hidden_layer_sizes=(50, 25), max_iter=1000, random_state=42)
# 使用整个数据集拟合神经网络分类器。
mlp.fit(data_transposed, numeric_labels)

# 计算基于权重的特征重要性。
# 对于每个输入特征（基因），计算与第一个隐藏层中的所有节点的权重的绝对值之和。
feature_importances = np.abs(mlp.coefs_[0]).sum(axis=1)
# 获取基因的索引，按特征重要性从高到低排序。
sorted_gene_indices = feature_importances.argsort()[::-1]
# 使用排序后的索引提取基因名。
sorted_genes = data_transposed.columns[sorted_gene_indices]
# 选择前 15 个最重要的基因。
selected_genes_nn = sorted_genes[:15].tolist()

# 为前 15 个基因绘制特征重要性条形图。
plt.figure(figsize=(10, 7))
# 绘制条形图，注意[::-1] 是为了确保最重要的基因显示在顶部。
plt.barh(selected_genes_nn[::-1], feature_importances[sorted_gene_indices][:15][::-1], color="skyblue")
plt.xlabel("Feature Importance")
plt.ylabel("Genes")
plt.title("Top 15 Genes Based on Neural Network Feature Importance")
plt.tight_layout()
# 将图像保存到指定路径。
plt.savefig("D:/BaiduNetdiskDownload/Diagnostic_video/136Diagnostic/13.SVM/U_D/py_nn_feature/NN_Feature_Importance.pdf")
plt.close()

# 将选定的基因写入一个新的文本文件。
with open("D:/BaiduNetdiskDownload/Diagnostic_video/136Diagnostic/13.SVM/U_D/py_nn_feature/NN_selected_genes.txt", 'w') as file:
    for gene in selected_genes_nn:
        file.write(f"{gene}\n")
