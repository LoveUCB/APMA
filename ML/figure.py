import warnings
warnings.filterwarnings('ignore')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from sklearn.preprocessing import LabelEncoder
encoder = LabelEncoder()
from sklearn import svm
import joblib
from xgboost import XGBClassifier
from lightgbm import LGBMClassifier
from catboost import CatBoostClassifier
from sklearn.linear_model import Lasso
from sklearn.linear_model import LassoCV
from sklearn.model_selection import StratifiedKFold
from sklearn.feature_selection import RFECV
from sklearn.feature_selection import RFE
from sklearn.model_selection import cross_val_predict
from sklearn.ensemble import VotingClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import AdaBoostClassifier
from sklearn.ensemble import StackingClassifier
from sklearn.metrics import accuracy_score
from sklearn.model_selection import cross_val_score
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import LabelEncoder
from tqdm import tqdm
from .model import model_rfe
from .model import stacking_model
from scipy.stats import wilcoxon

def plot_roc_curve(fpr, tpr,auc,fliename):
    """
    绘制ROC曲线
    参数：
    - fpr: 假正例率（False Positive Rate）列表
    - tpr: 真正例率（True Positive Rate）列表
    """
    plt.figure(figsize=(8,8))
    plt.plot(fpr, tpr, color='darkorange', lw=2, label='ROC curve')
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('(ROC) Curve')
    plt.legend(loc="lower right")
    plt.text(0.5, 0.3, 'AUC = {:.2f}'.format(auc), fontsize=12, ha='center')
    plt.savefig(fliename)
    plt.close()

def plot_spearman(input_file, output_folder):
    '''
    plot spearman correlation between the features
    :param input_file:
    :param output_folder:
    :return:
    '''
    with open(input_file, 'r') as f:
        columns = f.readline().strip().split('\t')
    data = pd.read_csv(input_file,skiprows=1,sep='\t')
    data = data.iloc[:,3:]
    data.columns = columns[3:]
    spearman_corr = data.corr(method="spearman")
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        print("Created")

    plt.figure(figsize=(12,10))
    sns.heatmap(spearman_corr,cmap="RdBu_r",fmt=".2f",
                annot=True,xticklabels='auto')
    plt.xticks(rotation=45,ha='right')
    plt.title("Spearman Correlation")

    output_path = os.path.join(output_folder, "spearman_corr.pdf")
    plt.savefig(output_path,format="pdf",bbox_inches='tight')
    plt.close()
    print(f"Spearman correlation plot generated")

def save_bar_chart_as_pdf(df,filename):
    """
    保存柱状图为PDF文件

    参数：
    categories: list，类别列表
    values: list，值列表
    filename: str，要保存的文件名
    """
    #print("saving bar chart... ",end = '')
    cores = [
        #svm.SVC(kernel="linear",max_iter=1000000),
        RandomForestClassifier(n_estimators=1000),
        GradientBoostingClassifier(n_estimators=1000),
        XGBClassifier(n_estimators=1000),
        LGBMClassifier(verbose=-1, n_estimators=1000)
    ]
    exp = [
        "Random Forest",
        "Gradient Boosting",
        "XGBoost",
        "LGBM"
    ]
    X = df.drop(columns=['Disease',"Site","Mutation"])
    y = encoder.fit_transform(df["Disease"])
    for i in range(len(cores)):
        cores[i].fit(X,y)
        values = cores[i].feature_importances_
        categories = X.columns
        sorted_data = sorted(zip(values, categories))
        values, categories = zip(*sorted_data)
        values = list(values)[::-1]
        categories = list(categories)[::-1]
        plt.figure(figsize=(12, 12))
        plt.bar(categories, values)
        plt.xticks(rotation=45, ha='right')
       #plt.xlabel('Categories')
        plt.ylabel('Values')
        plt.title('Feature Importances')
        plt.savefig(filename + "_" + exp[i] + ".pdf", format='pdf')
        plt.close()
    #print("Done")

def plot_roc_for_disease_pairs(file_path, output_dir):
    """
    Plot ROC curves for each pair of diseases with all features.

    Parameters:
    file_path (str): Path to the input data file.
    output_dir (str): Directory to save the output PDF files.

    Returns:
    None
    """
    # Read the txt file
    data = pd.read_csv(file_path, delimiter='\t')
    data = data.drop("Site", axis = 1)
    data = data.drop("Mutation", axis = 1)
    # Get unique disease categories
    diseases = data['Disease'].unique()

    # Generate pairs of diseases
    disease_pairs = [(diseases[i], diseases[j]) for i in range(len(diseases)) for j in range(i + 1, len(diseases))]
    # print(disease_pairs)
    # Plot for each disease pair
    for disease_pair in disease_pairs:
        # Create a figure and axis
        data_current = data[data["Disease"].isin(list(disease_pair))]
        fig, ax = plt.subplots(figsize=(10, 10))
        ax.plot([0, 1], [0, 1], color='grey', lw=1.5, linestyle='--')
        # Dictionary to store AUC values for each feature
        auc_dict = {}
        # Plot ROC curves for each feature
        for feature in data.columns[1:]:
            # Extract feature data and labels
            feature_data = data_current[[feature, 'Disease']].copy()
            disease_counts = feature_data['Disease'].value_counts()
            # 找到出现频次较多的疾病名称
            most_common_disease = disease_counts.idxmax()
            # 将出现频次较多的疾病名称设为0，其他疾病名称设为1
            feature_data['Disease'] = feature_data['Disease'].apply(lambda x: 1 if x == most_common_disease else 0)
            X = feature_data[[feature]]
            y = feature_data['Disease']
            # print(y)
            # Compute ROC curve
            fpr, tpr, _ = roc_curve(y, X)
            roc_auc = auc(fpr, tpr)
            # if auc <= 0.5 then reverse the y
            if roc_auc <= 0.5:
                y = [0 if m == 1 else 1 for m in y]
            fpr, tpr, _ = roc_curve(y, X)
            roc_auc = auc(fpr, tpr)
            # Plot ROC curve
            ax.plot(fpr, tpr, label=f'{feature} (AUC = {roc_auc:.4f})')
            # Store AUC value
            auc_dict[feature] = roc_auc

        # Sort legend labels by AUC values
        handles, labels = ax.get_legend_handles_labels()
        labels_and_aucs = [(label, auc_dict[label.split()[0]]) for label in labels]
        # print(labels_and_aucs)
        labels_and_aucs_sorted = sorted(labels_and_aucs, key=lambda x: x[1], reverse=True)
        labels_sorted = [x[0] for x in labels_and_aucs_sorted]
        handles_sorted = [handles[labels.index(label)] for label in labels_sorted]
        ax.legend(handles_sorted, labels_sorted, loc='lower right',fontsize='small')

        # Set title and axis labels
        plt.title(f'ROC Curves for {disease_pair[0]} vs {disease_pair[1]}')
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')

        # Save as PDF file
        output_path = f'{output_dir}/{disease_pair[0]}_vs_{disease_pair[1]}_ROC.pdf'
        plt.savefig(output_path, format='pdf')
        # Close the figure
        plt.close(fig)
    print("ROC plots generated")



def plot_box(data_file, output_folder):
    '''
    plot all features boxplot
    '''
    # 读取数据
    data = pd.read_csv(data_file,sep='\t')
    data = data.drop("Site", axis = 1)
    data = data.drop("Mutation", axis = 1)
    # 获取第一列疾病名称
    diseases = data.iloc[:, 0].unique()
    os.makedirs(output_folder, exist_ok=True)
    
    # 循环处理每个特征
    for col_index in range(1, len(data.columns)):
        feature_name = data.columns[col_index]
        plt.figure(figsize=(6, 6))
        
        for disease in diseases:
            # 获取特定疾病的数据
            disease_data = data[data.iloc[:, 0] == disease].iloc[:, col_index]
            # 生成水平坐标，抖散散点
            jittered_positions = np.random.normal(diseases.tolist().index(disease), 0.08, size=len(disease_data))
            # 绘制箱线图
            plt.boxplot(disease_data, positions=[diseases.tolist().index(disease)],showfliers=False,widths=0.4)
            # 绘制带透明度的散点图
            plt.scatter(jittered_positions, disease_data, alpha=0.5)
            plt.violinplot(disease_data,positions=[diseases.tolist().index(disease)],widths=0.6)
        
        
        plt.xticks(range(len(diseases)), diseases)
        plt.xlabel('Disease')
        # plt.ylabel(feature_name)
        plt.title(f'{feature_name}')
        plt.tight_layout()
        
        # 保存为PDF文件
        output_file = os.path.join(output_folder, f'{feature_name}.pdf')
        plt.savefig(output_file)
        plt.close()
    print("Boxplot generated")

def plot_importence_bar(df,filename):
    """
    保存柱状图为PDF文件

    参数：
    categories: list，类别列表
    values: list，值列表
    filename: str，要保存的文件名
    """
    df = pd.read_csv(df,sep='\t')
    print("saving bar chart... ",end = '')
    cores = [
        #svm.SVC(kernel="linear",max_iter=1000000),
        RandomForestClassifier(n_estimators=1000),
        #GradientBoostingClassifier(n_estimators=1000),
        XGBClassifier(n_estimators=1000),
        LGBMClassifier(verbose=-1, n_estimators=1000)
    ]
    exp = [
        "Random Forest",
        #"Gradient Boosting",
        "XGBoost",
        "LGBM"
    ]
    X = df.drop(columns=['Disease',"Site","Mutation"])
    y = encoder.fit_transform(df["Disease"])
    for i in range(len(cores)):
        cores[i].fit(X,y)
        values = cores[i].feature_importances_
        categories = X.columns
        sorted_data = sorted(zip(values, categories))
        values, categories = zip(*sorted_data)
        values = list(values)[::-1]
        categories = list(categories)[::-1]
        plt.figure(figsize=(9, 11))
        plt.bar(categories, values)
        plt.xticks(rotation=45, ha='right')
        plt.xlabel('Categories')
        plt.ylabel('Values')
        plt.title('Feature Importances')
        plt.savefig(filename + "/Importance_" + exp[i] + ".pdf", format='pdf')
        plt.close()
    print("Done")





