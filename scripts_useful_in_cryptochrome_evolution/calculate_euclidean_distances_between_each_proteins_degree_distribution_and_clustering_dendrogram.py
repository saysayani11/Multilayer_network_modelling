import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram
import os
os.chdir(r"C:\Users\saySa\OneDrive\Desktop\task_1\PSNs\networks\Degree")

def read_data(file_name):
    """Reads the .xlsx data into a pandas DataFrame."""
    return pd.read_excel(file_name, index_col=0)

def compute_pairwise_distance(df):
    """Compute pairwise Euclidean distances between columns of the DataFrame."""
    distances = pdist(df.T, metric='euclidean')
    distance_matrix = squareform(distances)
    return distance_matrix

def save_distance_data(distance_matrix, output_file):
    """Save the distance matrix as .xlsx."""
    df = pd.DataFrame(distance_matrix)
    df.to_excel(output_file)

def perform_hierarchical_clustering(distance_matrix):
    """Performs hierarchical clustering and returns linkage matrix."""
    # Note: 'ward' method minimizes the variance of the distances between clusters
    return linkage(distance_matrix, method='ward')

def plot_dendrogram(Z, labels, output_file):
    """Plot and save dendrogram."""
    plt.figure(figsize=(10, 7))
    dendrogram(Z, labels=labels, orientation='right')
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.show()

def main():
    input_file = 'Degree_Distribution.xlsx'
    output_distance_file = 'pairwise_distances.xlsx'
    output_dendrogram_file = 'dendrogram.png'
    
    # Read data
    data = read_data(input_file)
    
    # Compute pairwise Euclidean distances
    distance_matrix = compute_pairwise_distance(data)
    
    # Save the distance matrix
    save_distance_data(distance_matrix, output_distance_file)
    
    # Perform hierarchical clustering
    Z = perform_hierarchical_clustering(distance_matrix)
    
    # Plot and save dendrogram
    labels = data.columns
    plot_dendrogram(Z, labels, output_dendrogram_file)

if __name__ == '__main__':
    main()
