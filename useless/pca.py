# -*- coding: utf-8 -*-

"""

@ author: Jingran Wang

@ Email: jrwangspencer@stu.suda.edu.cn

@ Address: Center for Systems Biology, Department of Bioinformatics, School of Biology and Basic Medical Sciences, Soochow University, Suzhou 215123, China.

@ GitHub: https://github.com/Spencer-JRWang/APMA


"""

#############################################
### Introduction of PCA module
#
# @ This module is to excute PCA on the mutations
#
#############################################

import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import numpy as np
from matplotlib.patches import Ellipse
import matplotlib.colors as mcolors

class PCAPlotter:
    def __init__(self, data_file_path, output_file_path, color_list = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']):
        """
        Initializes PCAPlotter object.

        Parameters:
        - data_file_path (str): File path of the data.
        - output_file_path (str): File path to save the plot.
        """
        self.data_file_path = data_file_path
        self.output_file_path = output_file_path
        self.color_list = color_list
        self.data = None
        self.labels = None
        self.features = None
        self.principal_df = None

    def read_data(self):
        """
        Read data from the specified file.
        """
        # Read data
        self.data = pd.read_csv(self.data_file_path, sep='\t')
        # Drop 'Mutation' and 'Site' columns
        self.data = self.data.drop(['Mutation', 'Site'], axis=1)
        # Extract labels and features
        self.labels = self.data.iloc[:, 0]
        self.features = self.data.iloc[:, 1:]

    def standardize_data(self):
        """
        Standardize the features.
        """
        # Standardize the features
        scaler = StandardScaler()
        self.scaled_features = scaler.fit_transform(self.features)

    def perform_pca(self, n_components=2):
        """
        Perform Principal Component Analysis (PCA).

        Parameters:
        - n_components (int): Number of components for PCA.
        """
        # Perform PCA analysis
        pca = PCA(n_components=n_components)
        principal_components = pca.fit_transform(self.scaled_features)
        self.principal_df = pd.DataFrame(data=principal_components, columns=[f'PC{i+1}' for i in range(n_components)])
        # Add labels to the principal components
        self.principal_df['Label'] = self.labels

    def plot_ellipse(self, ax, data, label, color, n_std=2.0, **kwargs):
        """
        Plot covariance ellipses for specified data.

        Parameters:
        - ax (matplotlib.axes.Axes): Axes object to plot on.
        - data (numpy.ndarray): Data for which to plot ellipse.
        - label (str): Label for the data.
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

    def plot_pca(self):
        """
        Plot PCA scatter plot with covariance ellipses and mutation labels.
        """
        plt.figure(figsize=(7, 7))
        ax = plt.gca()
        if self.color_list:
            colors = self.color_list
        else:
            colors = plt.cm.get_cmap('tab10', len(self.principal_df['Label'].unique()))

        import itertools
        markers = itertools.cycle(['o', 's', '^', 'D', '*', 'x', '+', 'v', '<', '>'])

        for idx, label in enumerate(self.principal_df['Label'].unique()):
            subset = self.principal_df[self.principal_df['Label'] == label]
            color = colors[idx] if isinstance(colors, list) else colors(idx)
            marker = next(markers)
            plt.scatter(subset.iloc[:, 0], subset.iloc[:, 1], label=label, marker=marker, color = color, s=18)
            self.plot_ellipse(ax, subset.iloc[:, :2].values, label=label, color=color)



        plt.xlabel('PCA 1')
        plt.ylabel('PCA 2')
        plt.title('PCA of Mutations')
        plt.legend()
        plt.savefig(self.output_file_path + '/PCA.pdf', format='pdf')
        plt.close()

    def execute(self):
        """
        Execute the PCA analysis and plot generation.
        """
        self.read_data()
        self.standardize_data()
        self.perform_pca()
        self.plot_pca()

if __name__ == "__main__":
    # Usage
    data_file_path = '/Users/wangjingran/Downloads/APMA_outcome-28/paras.txt'
    output_file_path = '/Users/wangjingran/Desktop/SpencerW-APMA/APMA/Outcome'
    plotter = PCAPlotter(data_file_path, output_file_path)
    plotter.execute()
