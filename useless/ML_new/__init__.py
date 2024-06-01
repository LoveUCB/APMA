import warnings
warnings.filterwarnings('ignore')
import pandas as pd
import numpy as np
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

def ML_Build(category,file = '/Users/wangjingran/Desktop/SpencerW-APMA/APMA/data/paras.txt'):
    '''
    This function is intended to search for the best model and feature combination
    using machine learning algorithms based on given data file (default: paras.txt).
    Inputs: 
        - file : a string indicating the path of the input parameter file;
     Outputs:
         - The best feature combination and model combination found by searching are printed out in console.
           A txt file named "best_params.csv" will be generated with these information.
    '''
    print("...Loading data...")
    df_all = pd.read_csv(file, sep='\t')
    print("Success")
    cores = [
        #svm.SVC(kernel="linear",max_iter=1000000),
        #RandomForestClassifier(n_estimators=3000),
        #GradientBoostingClassifier(n_estimators=1000),
        #XGBClassifier(n_estimators = 1000)#,
        LGBMClassifier(verbose=-1, n_estimators=1000)
             ]

    exp = [
        #"SVM",
        #"RandomForest",
        #"GradientBoost",
        #"XGBoost"#,
        "LightGBM"
            ]




    from itertools import combinations
    from .explain import model_explain
    category = list(set(category))
    category = list(combinations(category, 2))
    from .figure import plot_roc_curve
    from .figure import save_bar_chart_as_pdf
    RFE_outcome = []
    f = open("/Users/wangjingran/Desktop/SpencerW-APMA/APMA/Outcome/Feature_selection.txt","w")
    # 如果category在搜索的字段中就进行标签归类
    for i in category:
        if i[0] == "Control" or i[0] == "control" or i[0] == "Control":
            Cat_A = i[0]
            Cat_B = i[1]

        elif i[0] == "Disease" or i[0] == "Severe" or i[0] == "disease" or i[0] == "severe":
            Cat_A = i[1]
            Cat_B = i[0]

        elif i[1] == "Control" or i[1] == "control":
            Cat_A = i[1]
            Cat_B = i[0]

        elif i[1] == "Disease" or i[1] == "Severe" or i[1] == "disease" or i[1] == "severe":
            Cat_A = i[0]
            Cat_B = i[1]

        else:
            Cat_A, Cat_B = i[0], i[1]

        f.write("--------------------" + str(Cat_A) + " and " + str(Cat_B) + "--------------------\n")
        print("--------------------" + str(Cat_A) + " and " + str(Cat_B) + "--------------------")
        df = df_all[df_all['Disease'].isin([Cat_A, Cat_B])]
        Site = df.reset_index(drop=True)[['Site']]
        Mutation = df.reset_index(drop=True)[['Mutation']]
        print("...Saving Importance Barplot...")
        save_bar_chart_as_pdf(df,'/Users/wangjingran/Desktop/SpencerW-APMA/APMA/Outcome/Figure/Importance/Importance_'+ str(Cat_A) + " vs " + str(Cat_B))
        print("Success")
        df = df.drop(columns = ["Site","Mutation"])
        RFE_Cat = []
        for j in range(len(cores)):
            print("...Searching for best feature combination...")
            f.write("#######" + exp[j] + "#######\n")
            print("#######" + exp[j] + "#######")
            ot = model_rfe(f, cores[j], df, Cat_A, Cat_B)
            print("Success")
            RFE_Cat.append(ot)
            print("...Exlaining The Model...")
            shap_v = model_explain(exp[j], df[ot[0]], df['Disease'].map({Cat_A: 0, Cat_B: 1}), f"{Cat_A} vs {Cat_B}")
            print("Success")
            print("...Building Stacking Model...")
            Best_IndegratedScore = stacking_model(f,Site,Mutation, df[ot[0]], df['Disease'].map({Cat_A: 0, Cat_B: 1}))
            print("Success")
            fpr, tpr, thresholds = roc_curve(Best_IndegratedScore.iloc[:, 0], Best_IndegratedScore.iloc[:, 2])
            roc_auc = auc(fpr, tpr)
            # 保存文件和图像
            print("...Saving ROC Curve...")
            plot_roc_curve(fpr,tpr,roc_auc,'/Users/wangjingran/Desktop/SpencerW-APMA/APMA/Outcome/Figure/ROC/ML/' + exp[j] +"_"+ str(Cat_A) + " vs " + str(Cat_B)+".pdf")
            print("Success")
            print("...Saving Files...")
            Best_IndegratedScore.iloc[:, 0] = Best_IndegratedScore.iloc[:, 0].map({0: Cat_A, 1: Cat_B})
            Best_IndegratedScore.to_csv('/Users/wangjingran/Desktop/SpencerW-APMA/APMA/Outcome/Score/' + str(exp[j]) + "_" + str(Cat_A) + " vs " + str(Cat_B) +'.txt',
                sep='\t', index=False, header=True)
            print("Success")
        RFE_outcome.append(RFE_Cat)
    f.close()


