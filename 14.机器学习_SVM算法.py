# 导入所需的库
import pandas as pd
import numpy as np
from sklearn.feature_selection import RFECV
from sklearn.svm import SVC
from sklearn.model_selection import StratifiedKFold
import matplotlib.pyplot as plt

# 定义函数，该函数使用SVM-RFE进行特征选择并绘制结果。
def perform_svm_rfe_updated(data_path, image_save_path, features_save_path):
    # 使用pandas的read_csv函数从指定的路径加载数据文件。这个文件是一个以tab分隔的文本文件，其中第一列作为索引。
    data = pd.read_csv(data_path, sep="\t", index_col=0)
    # 将数据转置，使样本为行，基因为列。
    data_transposed = data.transpose()

    # 从行名中提取样本类型（例如“con”）作为分类标签。
    labels = data_transposed.index.str.split("_").str[-1]
    # 将标签转换为数字，其中“con”标签为0，其他标签为1。
    numeric_labels = np.where(labels == "con", 0, 1)

    # 初始化线性支持向量机（SVM）分类器。
    svc_linear = SVC(kernel="linear", C=1)

    # 使用线性SVM初始化RFECV（带交叉验证的递归特征消除）。
    # step=0.05表示每次迭代时要移除的特征比例，cv=StratifiedKFold(5)表示使用分层5折交叉验证。
    min_features_to_select = 1
    rfecv_linear = RFECV(estimator=svc_linear, step=0.05, cv=StratifiedKFold(5), scoring="accuracy", min_features_to_select=min_features_to_select)
    # 使用数据和标签拟合RFECV对象。
    rfecv_linear.fit(data_transposed, numeric_labels)

    # 从cv_results_中提取均值测试得分。
    mean_scores = rfecv_linear.cv_results_['mean_test_score']

    # 绘制均值交叉验证得分图。
    plt.figure(figsize=(8, 6))
    plt.plot(range(1, len(mean_scores) + 1), mean_scores, color="darkgreen", label="Mean Accuracy")
    # 找到得分最高的特征数量。
    optimal_num_features = np.argmax(mean_scores) + 1
    # 画一条垂直线表示最佳特征数量。
    plt.axvline(x=optimal_num_features, linestyle="--", color="blue", label=f"Optimal number of features = {optimal_num_features}")
    plt.xlabel("Number of Features")
    plt.ylabel("Mean Cross-Validation Accuracy")
    plt.title("SVM-RFE Feature Selection")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    # 保存图像到指定的路径。
    plt.savefig(image_save_path)
    plt.close()

    # 从RFECV对象中提取选定的特征。
    selected_features = data_transposed.columns[rfecv_linear.support_].tolist()
    # 将选定的特征保存到一个新的文本文件中。
    with open(features_save_path, 'w') as file:
        for feature in selected_features:
            file.write(f"{feature}\n")

    return selected_features

# 以下是调用上述函数的代码。
# 替换以下路径为您的实际文件路径。
data_path = "D:/BaiduNetdiskDownload/Diagnostic_video/136Diagnostic/13.SVM/U_D/py_SVM/diffGeneExp.txt"
image_save_path = "D:/BaiduNetdiskDownload/Diagnostic_video/136Diagnostic/13.SVM/U_D/py_SVM/SVM-RFE.pdf"
features_save_path = "D:/BaiduNetdiskDownload/Diagnostic_video/136Diagnostic/13.SVM/U_D/py_SVM/SVM-RFE_selected_features.txt"

# 调用函数。
selected_features = perform_svm_rfe_updated(data_path, image_save_path, features_save_path)
