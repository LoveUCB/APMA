import warnings
warnings.filterwarnings('ignore')
import os
import shutil
import pandas as pd
import numpy as np
from sklearn import svm
from xgboost import XGBClassifier
from lightgbm import LGBMClassifier
from catboost import CatBoostClassifier
from sklearn.linear_model import Lasso
from sklearn.linear_model import LassoCV
from sklearn.model_selection import GridSearchCV
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
from itertools import combinations



def model_rfe(f,core,df,cat_A,cat_B):
    X = df.drop("Disease",axis = 1)
    y = df['Disease']
    X = X.reset_index(drop=True)
    y = y.reset_index(drop=True)
    shuffle_index = np.random.permutation(X.index)
    X = X.iloc[shuffle_index]
    y = y.iloc[shuffle_index]
    y_encode = y.map({cat_A: 0, cat_B: 1})
    outcome_feature = []
    outcome_score = []
    #print("Best Feature Combination Detecting...",end=' ')
    for i in tqdm(range(X.shape[1])):
        rfe = RFE(core, n_features_to_select=i + 1)
        rfe.fit(X, y_encode)
        selected_features = X.columns[rfe.support_]
        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
        scores = cross_val_score(core, X[selected_features], y_encode, cv=cv)
        selected_features = X.columns[rfe.support_]
        outcome_feature.append(selected_features)
        outcome_score.append(scores.mean())
    
    # 数据层面的最佳组合
    max_predict_data = max(outcome_score)
    best_features = list(outcome_feature[outcome_score.index(max_predict_data)])
    print("Best Features Combination Detected: " + str(best_features))
    print("Best Validation Score: " + str(max_predict_data))
    f.write("Best Features Combination Detected: " + str(best_features) + "\n")
    f.write("Best Validation Score: " + str(max_predict_data) + "\n")

    # 手动特征调优
    if "Effectiveness" not in best_features \
    and "Sensitivity" not in best_features \
    and "Stiffness" not in best_features \
    and "MSF" not in best_features \
    and "DFI" not in  best_features:
        best_features.append("Effectiveness")
    
    # 保证ddG在整合模型的考虑范围内
    if "ddG" not in best_features:
        best_features.append("ddG")
    
    # 当Degree在并且Eigenvector不在的时候,把Degree换成Eigenvector
    if "Degree" in best_features \
    and "Eigenvector" not in best_features:
        best_features.append("Eigenvector")
        best_features.remove("Degree")
    elif "Degree" in best_features \
    and "Eigenvector" in best_features:
        best_features.remove("Degree")
    
    # 保证氨基酸接触网络的Betweenness
    if "Betweenness" not in best_features \
    and "Closeness" not in best_features \
    and "Degree" not in best_features \
    and "Eigenvector" not in best_features \
    and "Clustering_coefficient" not in best_features:
        best_features.append("Betweenness")
    
    # 手动调优后的模型组合再次进行交叉验证
    scores_adj = cross_val_score(core, X[best_features], y_encode, cv=cv)
    scores_adj = scores_adj.mean()
    print("The Integrated Features will be: " + str(best_features))
    print("Best Validation Score For Integrated Features: " + str(scores_adj))
    f.write("The Integrated Features will be: " + str(best_features) + "\n")
    f.write("Best Validation Score For Integrated Features: " + str(scores_adj) + "\n")
    
    return best_features, scores_adj


def stacking_model(f,site,mutation,X,y_encode):
    scores_st = []
    X = X.reset_index(drop=True)
    y_encode = y_encode.reset_index(drop=True)
    stratified_kfold = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    shuffle_index = np.random.permutation(X.index)
    X = X.iloc[shuffle_index]
    y_encode = y_encode.iloc[shuffle_index]

    base_models = [
        ('RandomForest',RandomForestClassifier(n_estimators=2500, max_depth = 5)),
        ('GradientBoost',GradientBoostingClassifier(n_estimators=1000, max_depth = 5)),
        ('LGBM',LGBMClassifier(verbose = -1,n_estimators=1000, max_depth = 5)),
        ('XGBoost',XGBClassifier(n_estimators = 1000, max_depth = 5)),
        ('CatBoost',CatBoostClassifier(verbose = False,iterations = 1000, depth=5))
    ]
    meta_model = LogisticRegression(max_iter=10000000)

    # fetch probabilities
    print("...Training...")
    model_probabilities = {}
    for name, model in tqdm(base_models):
        probabilities = cross_val_predict(
                                            model, 
                                            X, 
                                            y_encode, 
                                            method='predict_proba', 
                                            cv = stratified_kfold
                                        )
        model_probabilities[name] = probabilities
    print("Success")
    
    AUCs = []
    Scores = []
    best_model_combination = []
    best_AUC = 0
    for r in range(1, len(base_models) + 1):
        for model_combination in combinations(base_models, r):
            # 按照组合形式构建特征矩阵
            features = np.concatenate([model_probabilities[name] for name, _ in model_combination], axis=1)
            meta_score = cross_val_predict(
                                            meta_model, 
                                            features, 
                                            y_encode, 
                                            method="predict_proba", 
                                            cv=stratified_kfold
                                           )
            Scores.append(meta_score)
            fpr, tpr, thresholds = roc_curve(y_encode, [i[1] for i in meta_score])
            roc_auc = auc(fpr, tpr)
            AUCs.append(roc_auc)
            print("Model: " + str([o[0] for o in model_combination]))
            print("AUC = " + str(roc_auc))
            f.write("Model: " + str([o[0] for o in model_combination]) + "\n")
            f.write("AUC = " + str(roc_auc) + "\n")
            # 找到最佳的模型组合
            if roc_auc > best_AUC:
                best_AUC = roc_auc
                best_model_combination = model_combination
    print("Best Stacking Model detected " + str(i[0] for i in best_model_combination))
    print("Best IntegratedScore AUC = " + str(max(AUCs)))
    f.write("Best Stacking Model detected " + str(i[0] for i in best_model_combination) + "\n")
    f.write("Best IntegratedScore AUC = " + str(max(AUCs)) + "\n")
    IndegratedScore = Scores[AUCs.index(max(AUCs))]
    scores_st = IndegratedScore[:, 1].tolist()
    dff = y_encode.to_frame()
    dff["Site"] = site.iloc[shuffle_index]
    dff["IntegratedScore"] = scores_st
    dff["Mutation"] = mutation.iloc[shuffle_index]
    return dff





