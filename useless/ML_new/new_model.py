from itertools import combinations
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from xgboost import XGBClassifier
from lightgbm import LGBMClassifier
from catboost import CatBoostClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_predict
from sklearn.datasets import make_classification
from sklearn.model_selection import cross_val_score
import numpy as np
from tqdm import tqdm

# 创建数据集
X, y = make_classification(n_samples=1000, n_features=20, random_state=42)

# 定义基础模型
base_models = [
        ('RandomForest',RandomForestClassifier(n_estimators=2500)),
        ('GradientBoost',GradientBoostingClassifier(n_estimators=1000)),
        ('LGBM',LGBMClassifier(verbose = -1,n_estimators=1000)),
        ('XGBoost',XGBClassifier(n_estimators = 1000)),
        ('CatBoost',CatBoostClassifier(verbose = False,iterations = 1000))
    ]

# 获取每个基础模型的预测概率
model_probabilities = {}
for name, model in base_models:
    probabilities = cross_val_predict(model, X, y, method='predict_proba', cv=5)
    model_probabilities[name] = probabilities

# 定义元模型
meta_model = LogisticRegression()

best_score = 0
best_model_combination = []

# 遍历所有模型组合
for r in range(1, len(base_models) + 1):
    for model_combination in combinations(base_models, r):
        # 按照组合形式构建特征矩阵
        features = np.concatenate([model_probabilities[name] for name, _ in model_combination], axis=1)
        # 交叉验证评估元模型
        scores = cross_val_score(meta_model, features, y, cv=5)
        mean_score = scores.mean()
        # 更新最佳得分和模型组合
        if mean_score > best_score:
            best_score = mean_score
            best_model_combination = model_combination

# 输出最佳模型组合和得分
print("Best model combination: ", [name for name, _ in best_model_combination])
print("Best score: ", best_score)
