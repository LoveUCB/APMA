# def accelerate_stacking_model(site, mutation, X, y_encode, meta_model = "Logestic Regression"):
#     if meta_model == "Logestic Regression":
#         meta = LogisticRegression(max_iter=10000000)
#     elif meta_model == "LASSO Regression":
#         meta = Lasso(alpha=0.01)
#     base_models = [
#         ('RandomForest',RandomForestClassifier(n_estimators=2500)),
#         ('GradientBoost',GradientBoostingClassifier(n_estimators=1000)),
#         ('LGBM',LGBMClassifier(verbose = -1,n_estimators=1000)),
#         ('XGBoost',XGBClassifier(n_estimators = 1000)),
#         ('CatBoost',CatBoostClassifier(verbose = False,iterations = 1000))
#     ]
#     stratified_kfold = StratifiedKFold(n_splits = 5, shuffle = True, random_state=42)
#     outcome = {
#         "RandomForest": [],
#         "GradientBoost": [],
#         "LGBM": [],
#         "XGBoost": [],
#         "CatBoost": []
#     }

#     # train the model
#     for i in base_models:
#         single_model = i[1]
#         single_model_exp = i[0]
#         score_model = cross_val_predict(single_model, X, y_encode, cv=stratified_kfold, method="predict_proba")
#         outcome[single_model_exp] = (score_model[:, 1])
    
#     # stacking model
#     all_model_exp = ["RandomForest", "GradientBoost", "LGBM", "XGBoost", "CatBoost"]
#     from itertools import combinations
#     all_combinations = []
#     for r in range(1, len(all_model_exp) + 1):
#         combinations_r = combinations(all_model_exp, r)
#         all_combinations.extend(combinations_r)

    
#     dff = y_encode.to_frame()
#     return dff
