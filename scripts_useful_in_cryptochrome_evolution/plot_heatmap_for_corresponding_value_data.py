import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt

# Set the working directory
os.chdir(r"C:\Users\saySa\OneDrive\Desktop\task_1\PDBeFOLD")


##################### PLOT HEATMAP FROM THE CLOSENESS DATA  #####################

# Set the style of the visualization
sns.set(style="white")

# Load the result DataFrame
result_df = pd.read_excel('result_dataset_betweenness.xlsx', engine='openpyxl')

# Convert the DataFrame to numeric, as it might contain non-numeric values like empty strings
result_df = result_df.apply(pd.to_numeric, errors='coerce')

# Create a heatmap from the DataFrame
plt.figure(figsize=(10, 8))
sns.heatmap(result_df, cmap="copper")

# Save the heatmap as a high-quality image
plt.savefig('heatmap.png', dpi=300, bbox_inches='tight')

# Display the heatmap
plt.show()
