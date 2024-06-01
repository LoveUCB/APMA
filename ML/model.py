# -*- coding: utf-8 -*-

"""

@ author: Jingran Wang

@ Email: jrwangspencer@stu.suda.edu.cn

@ Address: Center for Systems Biology, Department of Bioinformatics, School of Biology and Basic Medical Sciences, Soochow University, Suzhou 215123, China.

@ GitHub: https://github.com/Spencer-JRWang/APMA

"""

#############################################
### Introduction of model module
#
# @ This module is define ML models
#
#############################################



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
from sklearn.ensemble import StackingClassifier
from sklearn.metrics import accuracy_score
from sklearn.model_selection import cross_val_score
from tqdm import tqdm
from itertools import combinations

class ModelUtilities:
    def __init__(self):
        pass

    def model_rfe(self, f, core, df, cat_A, cat_B):
        """
        Perform recursive feature elimination (RFE) and return the best features and their validation score.

        Args:
        - f (file): File object for writing logs.
        - core: Classifier object for RFE.
        - df (DataFrame): DataFrame containing the dataset.
        - cat_A (str): Category A for encoding.
        - cat_B (str): Category B for encoding.

        Returns:
        - best_features (list): List of best features.
        - scores_adj (float): Best validation score for integrated features.
        """
        print("...RFE running...")
        X = df.drop("Disease", axis=1)
        y = df['Disease']
        X = X.reset_index(drop=True)
        y = y.reset_index(drop=True)
        shuffle_index = np.random.permutation(X.index)
        X = X.iloc[shuffle_index]
        y = y.iloc[shuffle_index]
        y_encode = y.map({cat_A: 0, cat_B: 1})
        outcome_feature = []
        outcome_score = []
        for i in tqdm(range(X.shape[1])):
            rfe = RFE(core, n_features_to_select=i + 1)
            rfe.fit(X, y_encode)
            selected_features = X.columns[rfe.support_]
            cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
            scores = cross_val_score(core, X[selected_features], y_encode, cv=cv)
            selected_features = X.columns[rfe.support_]
            outcome_feature.append(selected_features)
            outcome_score.append(scores.mean())

        max_predict_data = max(outcome_score)
        best_features = list(outcome_feature[outcome_score.index(max_predict_data)])
        f.write("Best Features Combination Detected: " + str(best_features) + "\n")
        f.write("Best Validation Score: " + str(max_predict_data) + "\n")
        print("Best Features Combination Detected: " + str(best_features))
        print("Best Validation Score: " + str(max_predict_data))

        if "Effectiveness" not in best_features \
        and "Sensitivity" not in best_features \
        and "Stiffness" not in best_features \
        and "MSF" not in best_features \
        and "DFI" not in  best_features:
            best_features.append("Effectiveness")

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

        if "Betweenness" not in best_features \
        and "Closeness" not in best_features \
        and "Degree" not in best_features \
        and "Eigenvector" not in best_features \
        and "Clustering_coefficient" not in best_features:
            best_features.append("Betweenness")

        scores_adj = cross_val_score(core, X[best_features], y_encode, cv=cv)
        scores_adj = scores_adj.mean()
        f.write("The Integrated Features will be: " + str(best_features) + "\n")
        f.write("Best Validation Score For Integrated Features: " + str(scores_adj) + "\n")
        print("The Integrated Features will be: " + str(best_features))
        print("Best Validation Score For Integrated Features: " + str(scores_adj))

        return best_features, scores_adj, outcome_score

    def model_combinations(self):
        """
        Generate all possible combinations of base models.

        Returns:
        - all_combinations (list of tuples): List containing all combinations of base models.
        """
        base_model = [
            #('RandomForest', RandomForestClassifier(n_estimators=2500)),
            ('GradientBoost', GradientBoostingClassifier(n_estimators=1000, max_depth=5)),
            ('LGBM', LGBMClassifier(verbose=-1, n_estimators=1000, max_depth=5)),
            ('XGBoost', XGBClassifier(n_estimators=1000, max_depth = 5)),
            ('CatBoost', CatBoostClassifier(verbose=False, iterations=1000, depth=5))
        ]
        all_combinations = []
        for r in range(1, len(base_model) + 1):
            combinations_r = combinations(base_model, r)
            all_combinations.extend(combinations_r)
        return all_combinations

    def stacking_model(self, site, mutation, X, y_encode, base_model):
        """
        Perform stacking ensemble and return the integrated scores.

        Args:
        - site (Series): Series containing site data.
        - mutation (Series): Series containing mutation data.
        - X (DataFrame): DataFrame containing feature data.
        - y_encode (Series): Series containing encoded target labels.
        - base_model (list of tuples): List of base models for stacking.

        Returns:
        - dff (DataFrame): DataFrame containing integrated scores.
        """
        scores_st = []
        X = X.reset_index(drop=True)
        y_encode = y_encode.reset_index(drop=True)
        stratified_kfold = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
        shuffle_index = np.random.permutation(X.index)
        X = X.iloc[shuffle_index]
        y_encode = y_encode.iloc[shuffle_index]
        meta_model = LogisticRegression(max_iter=10000000)
        stacking_clf = StackingClassifier(estimators=base_model, 
                                          final_estimator=meta_model, 
                                          stack_method='predict_proba'
                                        )
        score_st = cross_val_predict(stacking_clf, X, y_encode, cv=stratified_kfold, method="predict_proba")
        scores_st.append(score_st[:, 1])
        scores_st = np.array(scores_st)
        scores_st = np.mean(scores_st, axis=0)
        dff = y_encode.to_frame()
        dff["Site"] = site.iloc[shuffle_index]
        dff["IntegratedScore"] = scores_st
        dff["Mutation"] = mutation.iloc[shuffle_index]
        return dff


