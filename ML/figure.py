# -*- coding: utf-8 -*-

"""

@ author: Jingran Wang

@ Email: jrwangspencer@stu.suda.edu.cn

@ Address: Center for Systems Biology, Department of Bioinformatics, School of Biology and Basic Medical Sciences, Soochow University, Suzhou 215123, China.

@ GitHub: https://github.com/Spencer-JRWang/APMA

"""

#############################################
### Introduction of figure module
#
# @ This module is to draw basic figures
#
#############################################


import warnings
warnings.filterwarnings('ignore')
import pandas as pd
import numpy as np
import os
from sklearn.preprocessing import LabelEncoder
encoder = LabelEncoder()
from xgboost import XGBClassifier
from lightgbm import LGBMClassifier
from scipy.stats import ttest_ind
from scipy.stats import mannwhitneyu
import itertools
from sklearn.model_selection import cross_val_predict
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import GradientBoostingClassifier
from catboost import CatBoostClassifier
from sklearn.metrics import accuracy_score
from sklearn.model_selection import cross_val_score
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import LabelEncoder
from tqdm import tqdm
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from matplotlib.patches import Ellipse
import matplotlib.colors as mcolors

def plot_roc_curve(fpr, tpr, auc, filename):
    """
    Plot ROC curve with improved aesthetics and transparent area under the curve.

    Args:
        fpr (list): List of false positive rates.
        tpr (list): List of true positive rates.
        auc (float): Area under the ROC curve.
        filename (str): Name of the file to save the plot.
    """
    plt.figure(figsize=(5, 5))
    plt.plot(fpr, tpr, color='#1f77b4', lw=2, label='ROC curve (AUC = {:.2f})'.format(auc))
    plt.fill_between(fpr, tpr, color='#1f77b4', alpha=0.2)  # Transparent fill under the curve
    plt.plot([0, 1], [0, 1], color='gray', lw=1.5, linestyle='--', label='Random chance')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate', fontsize=10)
    plt.ylabel('True Positive Rate', fontsize=10)
    plt.title(filename.split('/')[-1])
    plt.legend(loc="lower right", fontsize=8)
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()


def plot_box(data_file, output_folder):
    '''
    Plot boxplots for all features and add t-test results between groups.

    Args:
        data_file (str): Path to the data file.
        output_folder (str): Path to the folder to save the boxplot PDF files.
    '''
    print("...Generating combined box plots...")
    # Read data
    data = pd.read_csv(data_file, sep='\t')
    data = data.drop(["Site", "Mutation"], axis=1)

    # Get unique disease names from the first column
    diseases = data.iloc[:, 0].unique()

    # Define the grid layout
    grid_rows = 4
    grid_cols = 4
    fig, axes = plt.subplots(grid_rows, grid_cols, figsize=(11, 11))
    # Loop through each feature
    for col_index, ax in enumerate(axes.flatten()):
        if col_index == 15:
            # Perform PCA and plot PCA scatter plot with ellipses
            plot_pca(data_file, ax)
            break

        feature_name = data.columns[col_index + 1]  # Skip the first column

        all_disease_data = []
        for disease in diseases:
            # Get data for a specific disease
            disease_data = data[data.iloc[:, 0] == disease].iloc[:, col_index + 1]
            all_disease_data.append(disease_data)

            # Generate jittered positions for scatter plot
            jittered_positions = np.random.normal(diseases.tolist().index(disease), 0.08, size=len(disease_data))
            # Plot boxplot
            ax.boxplot(disease_data, positions=[diseases.tolist().index(disease)], showfliers=False, widths=0.15,
                       patch_artist=True, boxprops=dict(facecolor="white"), medianprops=dict(color='black', linewidth=0.8))

            # Plot scatter plot with transparency
            ax.violinplot(disease_data, positions=[diseases.tolist().index(disease)], widths=0.35, showextrema=False,
                          bw_method='scott')
            ax.scatter(jittered_positions, disease_data, alpha=1, zorder=2, s=5, label=disease)

        # Perform t-tests between all disease pairs
        max_val = data.iloc[:, col_index + 1].max()
        spacing_factor = 0.2
        current_y = max_val + max_val * (spacing_factor-0.1)
        line_height = max_val * (spacing_factor-0.1) * 0.3  # Shorten the vertical line

        star_string = ''
        diff_1_list = []
        other_list = []
        for item in itertools.combinations(enumerate(diseases), 2):
            first_elements = [x[0] for x in item]
            if abs(first_elements[0] - first_elements[1]) == 1:
                diff_1_list.append(item)
            else:
                other_list.append(item)

        
        result = diff_1_list + other_list
        # obvious?
        obv_list = []
        for (i, disease1), (j, disease2) in result:
            if len(all_disease_data[i]) > 30 and len(all_disease_data[j]) > 30:
                # t
                t_stat, p_val = ttest_ind(all_disease_data[i], all_disease_data[j], equal_var=False)
            else:
                # Wilcoxon test
                _, p_val = mannwhitneyu(all_disease_data[i], all_disease_data[j])

            if p_val < 0.01:
                obv_list.append(1)
            else:
                obv_list.append(0)
        one_have = False
        
        for i in range(len(obv_list)):
            
            if i < len(diseases) - 1:
                
                if obv_list[i] == 1:
                    one_have = True
            
            else:
                pass

        if one_have:
            list_to_find_one = obv_list[:int(len(diseases) - 1)]
            list_to_find_one.reverse()
            index_of_one = list_to_find_one.index(1)
            last_index_of_one = len(list_to_find_one) - index_of_one - 1
            stop_index = result[last_index_of_one][1][0]

        for (i, disease1), (j, disease2) in result:
            if len(all_disease_data[i]) > 30 and len(all_disease_data[j]) > 30:
                # t
                t_stat, p_val = ttest_ind(all_disease_data[i], all_disease_data[j], equal_var=False)
            
            else:
                # Wilcoxon
                _, p_val = mannwhitneyu(all_disease_data[i], all_disease_data[j])

            if 0.001 <= p_val < 0.01:
                
                if j == i + 1:  # Adjacent groups
                    star_string = '**'
                    ax.plot([i, i, j, j], [current_y, current_y + line_height, current_y + line_height, current_y], lw=0.7, color='black')
                    ax.text((i + j) * 0.5, current_y, star_string, ha='center', va='bottom', fontsize=8, fontweight='bold')
                    
                    if j == stop_index:
                        current_y += line_height * 4
                
                else:  # Non-adjacent groups
                    star_string = '**'
                    ax.plot([i, i, j, j], [current_y, current_y + line_height, current_y + line_height, current_y], lw=0.7, color='black')
                    ax.text((i + j) * 0.5, current_y, star_string, ha='center', va='bottom', fontsize=8, fontweight='bold')
                    current_y += line_height * 4
            
            elif p_val < 0.001:
                
                if j == i + 1:  # Adjacent groups
                    star_string = '***'
                    ax.plot([i, i, j, j], [current_y, current_y + line_height, current_y + line_height, current_y], lw=0.7, color='black')
                    ax.text((i + j) * 0.5, current_y, star_string, ha='center', va='bottom', fontsize=8, fontweight='bold')
                    
                    if j == stop_index:
                        current_y += line_height * 4
                
                else:  # Non-adjacent groups
                    star_string = '***'
                    ax.plot([i, i, j, j], [current_y, current_y + line_height, current_y + line_height, current_y], lw=0.7, color='black')
                    ax.text((i + j) * 0.5, current_y, star_string, ha='center', va='bottom', fontsize=8, fontweight='bold')
                    current_y += line_height * 4

        # Customize plot appearance
        ax.set_xticks(range(len(diseases)))
        ax.set_xticklabels(diseases, rotation=45, ha='right')
        ax.set_ylabel(f'{feature_name}')
        
        if star_string:
            
            if col_index == 7:
                current_y += line_height * 4
                ax.set_ylim([None, current_y + max_val * (spacing_factor-0.05)])
            
            else:
                current_y -= line_height * 4
                ax.set_ylim([None, current_y + max_val * (spacing_factor-0.05)])
        
        else:
            pass

    # Adjust layout
    plt.tight_layout()

    # Save as a PDF file
    output_file = os.path.join(output_folder, 'combined_plots.pdf')
    plt.savefig(output_file)
    plt.close()


def plot_pca(data_file, ax):
    '''
    Perform PCA on the data and plot PCA scatter plot with ellipses.

    Args:
        data_file (str): Path to the data file.
        ax: Axis for plotting.
    '''
    data = pd.read_csv(data_file, sep='\t')
    data = data.drop(["Site", "Mutation"], axis=1)

    # Extract features
    features = data.iloc[:, 1:]

    # Standardize the features
    scaler = StandardScaler()
    scaled_features = scaler.fit_transform(features)

    # Perform PCA
    pca = PCA(n_components=2)
    principal_components = pca.fit_transform(scaled_features)
    principal_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'])

    # Add labels to the principal components
    principal_df['Label'] = data.iloc[:, 0]

    # Plot PCA scatter plot
    import itertools
    markers = itertools.cycle(['o', 's', '^', 'D', '*', 'x', '+', 'v', '<', '>'])
    colors = list(mcolors.TABLEAU_COLORS.values())  # Get a list of tableau colors
    for label, color in zip(principal_df['Label'].unique(), colors):
        subset = principal_df[principal_df['Label'] == label]
        marker = next(markers)
        ax.scatter(subset['PC1'], subset['PC2'], label=label, color=color, s = 5, marker = marker)

        # Plot ellipse for each category
        plot_ellipse(ax, subset[['PC1', 'PC2']].values, color)

    ax.set_xlabel('PC1', rotation=45, ha='right')
    ax.set_ylabel('PC2')
    # ax.grid(True, which='both', linestyle='--', linewidth=0.3, color='grey', alpha=0.3)
    ax.legend(fontsize = 5)

def plot_ellipse(ax, data, color, n_std=2.0, **kwargs):
    """
    Plot covariance ellipses for specified data.

    Parameters:
    - ax (matplotlib.axes.Axes): Axes object to plot on.
    - data (numpy.ndarray): Data for which to plot ellipse.
    - color (str): Color of the ellipse.
    - n_std (float): Number of standard deviations for ellipse size.
    """
    if data.size == 0:
        return
    # Calculate covariance and mean
    cov = np.cov(data, rowvar=False)
    mean = np.mean(data, axis=0)

    # Calculate eigenvalues and eigenvectors
    eigvals, eigvecs = np.linalg.eigh(cov)
    order = eigvals.argsort()[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]
    angle = np.degrees(np.arctan2(*eigvecs[:, 0][::-1]))

    # Calculate width and height of ellipse
    width, height = 2 * n_std * np.sqrt(eigvals)

    # Plot ellipse
    ell = Ellipse(xy=mean, width=width, height=height, angle=angle,
                  edgecolor=color, facecolor=mcolors.to_rgba(color, alpha=0.15), **kwargs)
    ax.add_patch(ell)
    return ell

def save_bar_chart_as_pdf(df, filename):
    """
    Save bar chart as a PDF file.

    Args:
        df (DataFrame): DataFrame containing the data.
        filename (str): Name of the file to save the PDF chart.
    """
    cores = [
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
    X = df.drop(columns=['Disease', "Site", "Mutation"])
    y = encoder.fit_transform(df["Disease"])
    for i in range(len(cores)):
        cores[i].fit(X, y)
        values = cores[i].feature_importances_
        categories = X.columns
        sorted_data = sorted(zip(values, categories))
        values, categories = zip(*sorted_data)
        values = list(values)[::-1]
        categories = list(categories)[::-1]
        plt.figure(figsize=(12, 12))
        plt.bar(categories, values)
        plt.xticks(rotation=45, ha='right')
        plt.ylabel('Values')
        plt.title('Feature Importances')
        plt.savefig(filename + "_" + exp[i] + ".pdf", format='pdf')
        plt.close()


def plot_roc_for_disease_pairs(file_path, output_dir):
    """
    Plot ROC curves for each pair of diseases with all features.

    Parameters:
    file_path (str): Path to the input data file.
    output_dir (str): Directory to save the output PDF files.

    Returns:
    None
    """
    print("...Generating feature roc plots...")
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
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.plot([0, 1], [0, 1], color='grey', lw=1.5, linestyle='--')
        # Dictionary to store AUC values for each feature
        auc_dict = {}

        # Plot ROC curves for each feature
        for feature in data.columns[1:]:

            # Extract feature data and labels
            feature_data = data_current[[feature, 'Disease']].copy()
            disease_counts = feature_data['Disease'].value_counts()

            # find most frequent disease
            most_common_disease = disease_counts.idxmax()

            # 0-1 encoding
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
        ax.legend(handles_sorted, labels_sorted, loc='lower right',fontsize=6)

        # Set title and axis labels
        plt.title(f'ROC Curves for {disease_pair[0]} vs {disease_pair[1]}')
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')

        # Save as PDF file
        output_path = f'{output_dir}/{disease_pair[0]}_vs_{disease_pair[1]}_ROC.pdf'
        plt.savefig(output_path, format='pdf')
        # Close the figure
        plt.close(fig)




def plot_importence_bar(df,filename):
    """
    Save bar chart of feature importances as PDF files.

    Args:
        df (str): Path to the DataFrame file.
        filename (str): Name of the folder to save the PDF files.
    """
    df = pd.read_csv(df,sep='\t')
    print("saving bar chart... ",end = '')
    cores = [
        #svm.SVC(kernel="linear",max_iter=1000000),
        RandomForestClassifier(n_estimators=1000),
        #GradientBoostingClassifier(n_estimators=1000),
        #XGBClassifier(n_estimators=1000),
        LGBMClassifier(verbose=-1, n_estimators=1000),
        CatBoostClassifier(verbose=False, iterations=1000)
    ]
    exp = [
        "Random Forest",
        #"Gradient Boosting",
        #"XGBoost",
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


def plot_spearman(input_file, output_folder):
    '''
    Plot Spearman correlation between the features.

    :param input_file: Path to the input file containing feature data.
    :param output_folder: Path to the folder where the output plot will be saved.
    :return: None
    '''
    # Read column names from the first line of the file
    with open(input_file, 'r') as f:
        columns = f.readline().strip().split('\t')
    
    # Read the data, skipping the first row (header)
    data = pd.read_csv(input_file, skiprows=1, sep='\t')
    data = data.iloc[:, 3:]  # Keep columns from the fourth one onwards
    data.columns = columns[3:]  # Assign the correct column names

    # Calculate Spearman correlation
    spearman_corr = data.corr(method="spearman")

    # Create the output folder if it does not exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        print("Created output folder:", output_folder)

    # Plotting the heatmap using Seaborn
    plt.figure(figsize=(12, 10))
    sns.heatmap(spearman_corr, cmap="RdBu_r", annot=True, fmt=".2f", annot_kws={"size": 10}, xticklabels='auto')
    plt.xticks(rotation=45, ha='right')
    plt.title("Feature Correlation")

    # Save the plot as a PDF file
    output_path = os.path.join(output_folder, "spearman.pdf")
    plt.savefig(output_path, format="pdf", bbox_inches='tight')
    plt.close()
    print(f"Spearman correlation plot saved to: {output_path}")


def plot_rfe(dict_data, output_folder):
    """
    Plots the RFE (Recursive Feature Elimination) results and saves the plot as a PDF file.

    Parameters:
    dict_data (dict): A dictionary where keys are labels (str) and values are lists of scores (list of float).
    output_folder (str): The directory where the output PDF file will be saved.

    """
    # Generate a list of counts from 1 to 15
    count = list(range(1, 16))
    
    # Set the figure size for the plot
    fig, ax = plt.subplots(figsize=(6, 6))

    # Define a list of markers to be used for different lines
    markers = ['o', 's', '^', 'D', '*', 'x', '+', 'v', '<', '>', 'o', 's', '^', 'D', '*', 'x', '+', 'v', '<', '>']
    i = 0  # Initialize the marker index
    
    # Define the fixed error value
    y_err = 0.008

    # Loop through each item in the dictionary
    for label, data in dict_data.items():
        # Plot the data with specified marker, line style, and width
        line, = ax.plot(count, data, marker=markers[i], linestyle='--', linewidth=1, label=label)

        # Define the fixed error value
        # y_err = [0.1 * value for value in data]
        
        # Plot the error bars at each point with the same color as the line
        ax.errorbar(count, data, yerr=y_err, fmt='none', ecolor=line.get_color(), capsize=2)
        
        # Find the index of the maximum value in the data
        max_index = data.index(max(data))
        
        # Highlight the maximum value with a purple marker
        ax.plot(count[max_index], data[max_index], 'purple', marker=markers[i])  # 'purple' with specified marker
        
        i += 1  # Increment the marker index

    plt.title('RFE Feature Selection')
    # Set the x-axis label
    plt.xlabel('Count')

    # Set the y-axis label
    plt.ylabel('Score')

    # Add a legend to the plot with specified font sizes
    plt.legend(fontsize=8, title_fontsize='8')

    # Save the plot as a PDF file in the specified output folder
    plt.savefig(f"{output_folder}/rfe.pdf", format='pdf')

    # Close the plot to free up memory
    plt.close()





def plot_dynamic_network(all_data_file, paras_file, output_file):
    """
    Plot parameter analysis based on input data files and save the plot as a PDF.

    Parameters:
    - all_data_file (str): File path to the main data file (e.g., '/path/to/alphafold_dyn.txt').
    - paras_file (str): File path to the parameters file (e.g., '/path/to/paras.txt').
    - output_file (str): File path to save the output PDF plot (e.g., '/path/to/parameters_analysis.pdf').
    """

    # Read data files
    all_data = pd.read_csv(all_data_file, sep="\t")
    paras = pd.read_csv(paras_file, sep="\t")

    # Rename columns for clarity
    all_data.rename(columns={all_data.columns[0]: "Site"}, inplace=True)

    # Extract unique residue and disease association
    site = paras[['Site', 'Disease']].drop_duplicates()

    # Define function to plot data for a specific parameter
    def plot_data(parameter, ax):
        """
        Plot data for a specific parameter.

        Parameters:
        - parameter (str): Parameter name to plot from 'all_data'.
        - ax (matplotlib.axes.Axes): Axes object to plot on.
        """
        # Extract data for the parameter
        df_parameter = all_data[['Site', parameter]]
        df_parameter.rename(columns={parameter: 'Value'}, inplace=True)
        
        # Merge with site data (residue and disease association)
        data = site.merge(df_parameter, on='Site')
        
        # Create base bar plot
        sns.barplot(x='Site', y='Value', data=df_parameter, color='#c5c5c5', edgecolor='none', width=1.2, ax=ax)
        
        # Add scatter plot with disease annotation
        sns.scatterplot(data=data, x='Site', y='Value', hue='Disease',
                        s=60, ax=ax)  # Use diamond marker
        
        # Set labels and titles
        ax.set_xlabel('Residues')
        ax.set_ylabel(parameter)
        residue_ticks = range(0, len(all_data) + 1, 100)
        ax.set_xticks(ticks=residue_ticks, labels=[str(i) for i in residue_ticks])
        ax.legend().remove()  # Remove legend

    # Parameters to plot
    parameters = ['Effectiveness', 'Sensitivity', 'Stiffness', 'DFI', 'MSF']

    # Create subplots for each parameter
    fig, axs = plt.subplots(nrows=len(parameters), figsize=(13, 13), sharex=True)

    # Plot each parameter
    for i, param in enumerate(parameters):
        plot_data(param, axs[i])

    # Adjust layout to avoid overlap
    plt.tight_layout()
    plt.legend()

    # Save plot as PDF
    plt.savefig(output_file, format='pdf')
    plt.close()


