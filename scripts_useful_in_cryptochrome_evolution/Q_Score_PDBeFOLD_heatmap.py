import os, pandas as pd, seaborn as sns, matplotlib.pyplot as plt
import numpy as np

os.chdir(r"C:\Users\saySa\OneDrive\Desktop\task_1\PDBeFOLD\RMSD")
df = pd.read_excel('rmsd_pdbefold.xlsx')

row_labels = [x.replace('.pdb:A', '') for x in df.iloc[1:, 1]]
col_labels = [x.replace('.pdb:A', '') for x in df.iloc[0, 2:]]

df_mod = df.iloc[1:, 2:].fillna(1).applymap(lambda x: pd.to_numeric(str(x).strip()))
df_mod.columns, df_mod.index = col_labels, row_labels

# Zero out the diagonal
np.fill_diagonal(df_mod.values, np.nan)

sns.set_style("white")
plt.figure(figsize=(10, 8), dpi=300)
sns.heatmap(df_mod, annot=True, cmap='RdBu', cbar_kws={"label": "RMSD"}, fmt='.2f', mask=np.isnan(df_mod))
plt.tight_layout()
plt.savefig("heatmap.png", dpi=300)
plt.show()
