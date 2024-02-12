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
    data = data.iloc[:,2:]
    data.columns = columns[2:]
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
    print(f"spearman correlation is saved to {output_path}")

def save_bar_chart_as_pdf(df,filename):
    """
    保存柱状图为PDF文件

    参数：
    categories: list，类别列表
    values: list，值列表
    filename: str，要保存的文件名
    """
    print("saving bar chart... ",end = '')
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
    X = df.drop(columns=['Disease',"Site"])
    y = encoder.fit_transform(df["Disease"])
    for i in range(len(cores)):
        cores[i].fit(X,y)
        values = cores[i].feature_importances_
        categories = X.columns
        sorted_data = sorted(zip(values, categories))
        values, categories = zip(*sorted_data)
        values = list(values)[::-1]
        categories = list(categories)[::-1]
        plt.figure(figsize=(16, 11))
        plt.bar(categories, values)
        plt.xticks(rotation=45, ha='right')
        plt.xlabel('Categories')
        plt.ylabel('Values')
        plt.title('Feature Importances')
        plt.savefig(filename + "_" + exp[i] + ".pdf", format='pdf')
        plt.close()
    print("Done")


