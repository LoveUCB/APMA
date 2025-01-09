# -*- coding: utf-8 -*-

"""

@ author: Jingran Wang

@ Email: jrwangspencer@stu.suda.edu.cn

@ Address: Center for Systems Biology, Department of Bioinformatics, School of Biology and Basic Medical Sciences, Soochow University, Suzhou 215123, China.

@ GitHub: https://github.com/Spencer-JRWang/APMA

"""

#############################################
### Introduction of ML module
#
# @ This module is to excute ML on the mutation features
#
#############################################



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
from .model import ModelUtilities
from .explain import model_explain

def ML_Build(category, file='/home/wangjingran/APMA/data/paras.txt'):
    '''
    This function searches for the best model and feature combination
    using machine learning algorithms based on the given data file (default: paras.txt).
    
    Inputs: 
        - file : a string indicating the path of the input parameter file;
        
    Outputs:
        - The best feature combination and model combination found by searching are printed out in the console.
          A txt file named "Model_connstruction.txt" will be generated with this information.
    '''
    
    # Read data from the given file
    print("...Machine learning starting...")
    df_all = pd.read_csv(file, sep='\t')
    
    # Define a list of ML models to be evaluated
    cores = [
        #svm.SVC(kernel="linear",max_iter=1000000),
        #RandomForestClassifier(n_estimators=3000),
        #GradientBoostingClassifier(n_estimators=1000),
        #XGBClassifier(n_estimators = 1000)#,
        LGBMClassifier(verbose=-1, n_estimators=1000, max_depth=5)
    ]

    # Define explanations for models
    exp = [
        #"SVM",
        #"RandomForest",
        #"GradientBoost",
        #"XGBoost"#,
        "LightGBM"
    ]

    # Import necessary modules
    from itertools import combinations
    from .explain import model_explain
    category = list(dict.fromkeys(category))
    category = list(combinations(category, 2))
    from .figure import plot_roc_curve
    # from .figure import save_bar_chart_as_pdf
    
    
    # Open a file to save feature selection outcome
    f = open("/home/wangjingran/APMA/Outcome/Model_construction.txt", "w")
    
    # Iterate through each combination of categories
    RFE_outcome = {}
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

        # Write category information to the feature selection file
        f.write("--------------------" + str(Cat_A) + " and " + str(Cat_B) + "--------------------\n")
        print("--------------------" + str(Cat_A) + " and " + str(Cat_B) + "--------------------")
        
        # Filter data based on categories
        df = df_all[df_all['Disease'].isin([Cat_A, Cat_B])]
        Site = df.reset_index(drop=True)[['Site']]
        Mutation = df.reset_index(drop=True)[['Mutation']]
        
        # Save bar chart as PDF
        # print("...Generating importance bar chart...")
        # save_bar_chart_as_pdf(df, '/Users/wangjingran/Desktop/APMA_dev/Outcome/Figure/Importance/Importance_' + str(Cat_A) + " vs " + str(Cat_B))
        df = df.drop(columns=["Site", "Mutation"])
        model_util = ModelUtilities()
        
        # Iterate through each model
        for j in range(len(cores)):
            f.write("#######" + exp[j] + "#######\n")
            print("#######" + exp[j] + "#######")
            ot = model_util.model_rfe(f, cores[j], df, Cat_A, Cat_B)
            RFE_outcome[f"{Cat_A} vs {Cat_B}"] = ot[2]
            print("...Model explaining...")
            shap_v = model_explain(exp[j], df[ot[0]], df['Disease'].map({Cat_A: 0, Cat_B: 1}), f"{Cat_A} vs {Cat_B}")
            all_com = model_util.model_combinations()
            AUCs = []
            Scores = []
            print("...Stacking model is building...")
            
            # Build stacking model
            for m in all_com:
                IntegratedScore = model_util.stacking_model(Site, Mutation, df[ot[0]], df['Disease'].map({Cat_A: 0, Cat_B: 1}), list(m))
                Scores.append(IntegratedScore)
                fpr, tpr, thresholds = roc_curve(IntegratedScore.iloc[:, 0], IntegratedScore.iloc[:, 2])
                roc_auc = auc(fpr, tpr)
                AUCs.append(roc_auc)
                f.write("Model: " + str([o[0] for o in m]) + "\n")
                f.write("AUC = " + str(roc_auc) + "\n")
                print("Model: " + str([o[0] for o in m]))
                print("AUC = " + str(roc_auc))
            print("====== Done ======")
            best_stacking = []
            for t in all_com[AUCs.index(max(AUCs))]:
                best_stacking.append(t[0])
            f.write("Best Stacking Model detected " + str(best_stacking) + "\n")
            f.write("Best IntegratedScore AUC = " + str(max(AUCs)) + "\n")
            print("Best Stacking Model detected " + str(best_stacking))
            print("Best IntegratedScore AUC = " + str(max(AUCs)))

            print("...Generating roc plots...")
            Best_IndegratedScore = Scores[AUCs.index(max(AUCs))]
            fpr, tpr, thresholds = roc_curve(Best_IndegratedScore.iloc[:, 0], Best_IndegratedScore.iloc[:, 2])
            roc_auc = auc(fpr, tpr)
            # Save files and images
            plot_roc_curve(fpr, tpr, roc_auc, '/home/wangjingran/APMA/Outcome/Figure/ROC/ML/' + exp[j] + "_" + str(Cat_A) + " vs " + str(Cat_B) + ".pdf")
            print("...Saving files...")
            Best_IndegratedScore.iloc[:, 0] = Best_IndegratedScore.iloc[:, 0].map({0: Cat_A, 1: Cat_B})
            Best_IndegratedScore.to_csv('/home/wangjingran/APMA/Outcome/Score/' + str(exp[j]) + "_" + str(Cat_A) + " vs " + str(Cat_B) + '.txt', sep='\t', index=False, header=True)

    f.close()

    print("...Generating rfe figures...")
    from .figure import plot_rfe
    plot_rfe(RFE_outcome, "/home/wangjingran/APMA/Outcome/Figure")



