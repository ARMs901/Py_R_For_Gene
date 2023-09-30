# 导入必要的库
import matplotlib.pyplot as plt
from matplotlib_venn import venn3

# 定义一个函数，用于从给定的文件路径中读取基因
def read_genes_from_file(file_path):
    # 使用'with'语句打开文件，这样在操作完成后文件会自动关闭
    with open(file_path, 'r') as f:
        # 读取文件的每一行，将其转换为一个集合并返回
        # 使用set是因为集合中的元素是唯一的，这对于基因列表来说是有用的
        return set([line.strip() for line in f.readlines()])

# 使用上面定义的函数读取提供的文件中的基因列表
genes_svm_rfe = read_genes_from_file("D:/BaiduNetdiskDownload/Diagnostic_video/136Diagnostic/14.venn/U_D/SVM-RFE.gene.txt")
genes_lasso = read_genes_from_file("D:/BaiduNetdiskDownload/Diagnostic_video/136Diagnostic/14.venn/U_D/LASSO.gene.txt")
genes_nn = read_genes_from_file("D:/BaiduNetdiskDownload/Diagnostic_video/136Diagnostic/14.venn/U_D/NN_selected_genes.txt")
#genes_svm_rfe = read_genes_from_file("D:/BaiduNetdiskDownload/Diagnostic_video/136Diagnostic/14.venn/U_D/SVM-RFE_selected_features.txt")

# 创建Venn图
plt.figure(figsize=(10, 7))  # 设置图的大小
# 使用venn3函数创建三个集合的Venn图，输入是三个基因集合，set_labels参数设置每个集合的标签
v = venn3([genes_svm_rfe, genes_lasso, genes_nn], set_labels=['SVM-RFE', 'LASSO', 'NN'])

# 计算三个基因列表的交集
intersected_genes = genes_svm_rfe.intersection(genes_lasso).intersection(genes_nn)
# 将交集基因写入文件
with open("D:/BaiduNetdiskDownload/Diagnostic_video/136Diagnostic/14.venn/U_D/intersected_genes.txt", 'w') as f:
    for gene in intersected_genes:
        f.write(f"{gene}\n")

# 保存Venn图为PNG格式到指定的位置
plt.savefig("D:/BaiduNetdiskDownload/Diagnostic_video/136Diagnostic/14.venn/U_D/venn_diagram.png")
# 保存Venn图为PDF格式到指定的位置
plt.savefig("D:/BaiduNetdiskDownload/Diagnostic_video/136Diagnostic/14.venn/U_D/venn_diagram.pdf")
# 在屏幕上显示Venn图
plt.show()
