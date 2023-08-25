import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

os.chdir(r"C:\Users\saySa\OneDrive\Desktop\task_1\PDBeFOLD")

df = pd.read_excel('rmsd_pdbefold.xlsx')

# Extracting the PDB IDs for rows and columns
row_labels = df.iloc[1:, 1].values
row_labels = [label.replace('.pdb:A', '') for label in row_labels]

col_labels = df.iloc[0, 2:].values
col_labels = [label.replace('.pdb:A', '') for label in col_labels]

df_modified = df.iloc[1:, 2:]
df_modified.fillna(1, inplace=True)
df_modified.columns = col_labels
df_modified.index = row_labels

# Strip whitespace and convert to numeric
df_modified = df_modified.applymap(lambda x: pd.to_numeric(str(x).strip(), errors='coerce'))

sns.set_style("white")
plt.figure(figsize=(10, 8), dpi=300)  # Set DPI here
sns.heatmap(df_modified, annot=True, cmap='RdBu', cbar_kws={"label": "RMSD Value"}, fmt='.2f')
plt.tight_layout()
plt.savefig("heatmap.png", dpi=300)

plt.show()
