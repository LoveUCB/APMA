from setuptools import setup, find_packages
setup(
    name='APMA',
    version='1.0',
    description='Auto Protein Mutation Analyzer',
    author='Wang Jingran',
    author_email='jrwangspencer@stu.suda.edu.cn',
    packages=find_packages(),
    include_package_data=True,
    url='https://github.com/LoveUCB/APMA',
    install_requires=[
        "pandas",
        "matplotlib",
        "numpy",
        "rpy2",
        "joblib",
        "prody",
        "seaborn",
        "catboost",
        "scikit-learn",
        "lightgbm",
        "xgboost",
        "Bio",
        "Biopython"
    ]
)
