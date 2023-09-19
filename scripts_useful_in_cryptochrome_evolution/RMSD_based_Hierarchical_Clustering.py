import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage

# Data preparation
os.chdir(r"C:\Users\saySa\OneDrive\Desktop\task_1\PDBeFOLD\RMSD")
df = pd.read_excel('RMSD_pdbefold.xlsx')

row_labels = [x.replace('.pdb:A', '') for x in df.iloc[1:, 1]]

df_mod = df.iloc[1:, 2:].fillna(1).applymap(lambda x: pd.to_numeric(str(x).strip()))
df_mod.index = row_labels

# Hierarchical clustering
linked = linkage(df_mod.fillna(0), 'single')

# Plot dendrogram
plt.figure(figsize=(14, 10), dpi=300)
dendrogram(linked, labels=df_mod.index, orientation='top')
plt.title('Hierarchical Clustering Dendrogram')
plt.xlabel('Cluster Labels')
plt.ylabel('RMSD Distance')
plt.show()
