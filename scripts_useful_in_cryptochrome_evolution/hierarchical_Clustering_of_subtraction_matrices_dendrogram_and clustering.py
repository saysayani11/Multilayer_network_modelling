import os
import numpy as np
import pandas as pd
from sklearn.cluster import AgglomerativeClustering, KMeans
from scipy.cluster.hierarchy import dendrogram, linkage
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import fcluster, dendrogram, linkage
from sklearn.metrics import silhouette_score

directory = r"C:\Users\saySa\OneDrive\Desktop\task_1\Distance_matrix_Analysis\matrix_substraction\data"

# Load matrices and flatten them, and store filenames
matrices_list = []
filenames_list = []  # List to store filenames
max_size = 0

for filename in os.listdir(directory):
    if filename.endswith(".csv"):
        file_path = os.path.join(directory, filename)
        matrix = pd.read_csv(file_path).values.flatten()
        matrices_list.append(matrix)
        filenames_list.append(filename)  # Save the filename
        max_size = max(max_size, matrix.size)

# Pad matrices to have the same number of elements
matrices_list = [np.pad(matrix, (0, max_size - matrix.size)) for matrix in matrices_list]

# Check if matrices_list is not empty
if not matrices_list:
    print("No matrices loaded. Please check the directory and file format.")
else:
    matrices_array = np.vstack(matrices_list)  # Convert list of 1D arrays to a 2D array
    
    # Hierarchical Clustering
    Z = linkage(matrices_array, 'ward')
    
    # Increase figure size and Plot Dendrogram
    plt.figure(figsize=(20, 10))  # Increase figure width
    plt.title('Hierarchical Clustering Dendrogram')
    dendrogram(Z)
    plt.xticks(rotation=45, ha='right', fontsize=5)  # Reduce font size
    plt.tight_layout()  # Adjust layout
    plt.savefig("dendrogram.png", dpi=300)
    plt.show()
    
############  COUNT THE NUMBER OF CLUSTERS BASED ON A CUTOFF VALUE ##############


# Assuming `matrices_array` is your data and `Z` is the linkage matrix
# Z = linkage(matrices_array, 'ward')

silhouette_scores = []

# Loop over some range of cluster numbers to compute silhouette scores
for num_clusters in range(2, 10):  # for example, checking from 2 to 9 clusters
    clusters = fcluster(Z, t=num_clusters, criterion='maxclust')
    score = silhouette_score(matrices_array, clusters)
    silhouette_scores.append(score)

# Find the optimal number of clusters
optimal_num_clusters = np.argmax(silhouette_scores) + 2  # Adding 2 as our range starts from 2

# Assign each matrix to its respective cluster
clusters = fcluster(Z, t=optimal_num_clusters, criterion='maxclust')

# Group matrices by cluster along with their filenames
clustered_matrices = {}
for i, (cluster_num, filename) in enumerate(zip(clusters, filenames_list)):
    if cluster_num not in clustered_matrices:
        clustered_matrices[cluster_num] = []
    clustered_matrices[cluster_num].append(f"Matrix {i+1} - {filename}")

# Visualize silhouette scores
plt.figure(figsize=(10, 6), dpi=300)
plt.plot(range(2, 10), silhouette_scores, marker='o')
plt.title('Silhouette Scores vs Number of Clusters')
plt.xlabel('Number of Clusters')
plt.ylabel('Silhouette Score')
plt.show()

print(f"The optimal number of clusters is: {optimal_num_clusters}")

with pd.ExcelWriter('clustered_matrices.xlsx') as writer:
    
    # Loop through each cluster and save the matrix data to a separate sheet
    for cluster_num, matrices in clustered_matrices.items():
        # Create a DataFrame for the current cluster
        df = pd.DataFrame(matrices, columns=['Matrix and Filename'])
        
        # Write the DataFrame to a new sheet in the Excel file
        df.to_excel(writer, sheet_name=f'Cluster {cluster_num}', index=False)
        
print("Clustered matrices data has been saved to 'clustered_matrices.xlsx'.")

