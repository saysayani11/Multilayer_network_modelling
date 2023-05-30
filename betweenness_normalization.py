import os
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

###--- Set Path
os.chdir(r"C:\Users\saySa\OneDrive\Desktop\test\betweenness")

correspondence_map = pd.read_excel('correspondence_map.xlsx')

###--- (a) Divide each column by its maximum value
correspondence_map = correspondence_map.apply(lambda x: x / x.max())

###--- (b) Drop the last column
correspondence_map.drop(correspondence_map.columns[-1], axis=1, inplace=True)

mask = correspondence_map.isna()

num_rows, num_cols = correspondence_map.shape
cell_width = 1
cell_height = 1
fig_width = cell_width * num_cols
fig_height = cell_height * num_rows
scale_factor = 0.1
fig_width *= scale_factor
fig_height *= scale_factor


plt.figure(figsize=(fig_width, fig_height), dpi=300)
cmap = "rocket"
sns.heatmap(correspondence_map, cmap=cmap, linewidths=0.5, linecolor='white', square=True)
protein_labels = ['3ZXS', '3UMV', '4GU5', '4U63', '2JDA']
plt.xticks(ticks=range(num_cols), labels=protein_labels, fontsize=4)
plt.show()
