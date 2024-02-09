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

def model_rfe(core,df,cat_A,cat_B):
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
    print("Best Feature Combination Detecting...",end=' ')
    for i in range(X.shape[1]):
        rfe = RFE(core, n_features_to_select=i + 1)
        rfe.fit(X, y_encode)
        selected_features = X.columns[rfe.support_]
        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
        scores = cross_val_score(core, X[selected_features], y_encode, cv=cv)
        selected_features = X.columns[rfe.support_]
        outcome_feature.append(selected_features)
        outcome_score.append(scores.mean())
    max_predict = max(outcome_score)
    print("Done")
    print("Best Features Combination Detected: " + str(list(outcome_feature[outcome_score.index(max_predict)])))
    print("Best Validation Score: " + str(max_predict))
    return list(outcome_feature[outcome_score.index(max_predict)]),max_predict

'''
Here maybe adding LASSOCV function
'''

'''
adding grid search
'''

def grid_search(X,y_encode):
    base_model = [
        ('RandomForest',RandomForestClassifier(n_estimators=2000)),
        ('GradientBoost',GradientBoostingClassifier(n_estimators=1500)),
        ('LGBM',LGBMClassifier(verbose = -1,n_estimators=1500)),
        ('XGBoost',XGBClassifier(n_estimators = 1500)),
        ('CatBoost',CatBoostClassifier(verbose = False,iterations = 800))
    ]
    from itertools import combinations
    all_combinations = []
    for r in range(1, len(base_model) + 1):
        combinations_r = combinations(base_model, r)
        all_combinations.extend(combinations_r)
    return all_combinations

def stacking_model(X,y_encode,base_model):
    scores_st = []
    X = X.reset_index(drop=True)
    y_encode = y_encode.reset_index(drop=True)
    stratified_kfold = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    shuffle_index = np.random.permutation(X.index)
    X = X.iloc[shuffle_index]
    y_encode = y_encode.iloc[shuffle_index]
    meta_model = meta_model = LogisticRegression(max_iter=10000000)
    stacking_clf = StackingClassifier(estimators=base_model, final_estimator=meta_model, stack_method='predict_proba')
    score_st = cross_val_predict(stacking_clf, X, y_encode, cv=stratified_kfold, method="predict_proba")
    scores_st.append(score_st[:, 1])
    scores_st = np.array(scores_st)
    scores_st = np.mean(scores_st, axis=0)
    dff = y_encode.to_frame()
    dff["IntegratedScore"] = scores_st
    return dff

def new_folder():
    grandparent_directory = os.path.abspath(os.path.join(os.getcwd(), os.pardir, os.pardir))
    new_folder_name = "IntegratedScore"
    new_folder_path = os.path.join(grandparent_directory, new_folder_name)
    if os.path.exists(new_folder_path):
        print("Old IntegratedScore will be deleted")
        shutil.rmtree(new_folder_path)
    new_folder_path = os.path.join(grandparent_directory, new_folder_name)
    os.makedirs(new_folder_path)

