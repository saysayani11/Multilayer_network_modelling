import os
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#--- Set Path
os.chdir(r"C:\Users\saySa\OneDrive\Desktop\test\betweenness")

correspondence_map = pd.read_excel('correspondence_map.xlsx')
mask = correspondence_map.isna()



# # Calculate mean betweenness centrality for each residue
residue_mean_centrality = correspondence_map.mean(axis=1)
# # residue_mean_centrality.to_excel("test.xlsx")

# # residues with high centrality above a threshold (e.g., top 10%)
# threshold = residue_mean_centrality.quantile(0.9)  # Example threshold of top 10%
# high_centrality_residues = residue_mean_centrality[residue_mean_centrality > threshold]

k1 = correspondence_map["Protein_1"]
k2 = residue_mean_centrality


df = pd.DataFrame({
    'Protein_1': k1,
    'Residue_Mean_Centrality': k2
})

df.dropna(subset = ['Protein_1'], inplace = True)



# threshold_new = df["Residue_Mean_Centrality"].quantile(0.9)
# high_centrality_residues2 = df[df["Residue_Mean_Centrality"] > threshold]