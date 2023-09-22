import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Specify the directory containing the .csv files
directory = r"C:\Users\saySa\OneDrive\Desktop\task_1\Distance_matrix_Analysis\matrix_substraction\data"

# Output directory to save the heatmaps
output_dir = r"C:\Users\saySa\OneDrive\Desktop\task_1\Distance_matrix_Analysis\matrix_substraction\single_heatmaps"
os.makedirs(output_dir, exist_ok=True)

# List all files in the directory
for filename in os.listdir(directory):
    # Check if the file is a .csv file
    if filename.endswith(".csv"):
        try:
            # Construct the full file path
            file_path = os.path.join(directory, filename)
            
            # Read the .csv file into a DataFrame
            df = pd.read_csv(file_path)
            
            # Set up the matplotlib figure
            plt.figure(figsize=(10, 8))
            
            # Get the number of rows and columns in the DataFrame
            nrows, ncols = df.shape
            
            # Decide interval for both x and y-axis labels for readability
            label_interval = max(1, min(ncols, nrows) // 10)
            
            # Set ticks at intervals
            x_ticks = np.arange(0, ncols, label_interval)
            y_ticks = np.arange(0, nrows, label_interval)
            
            # Draw the heatmap
            ax = sns.heatmap(df, cmap="viridis", cbar_kws={"shrink": .75})
            
            # Set and rotate the labels for better readability
            ax.set_xticks(x_ticks)
            ax.set_xticklabels(x_ticks + 1, rotation=45, horizontalalignment='right')
            ax.set_yticks(y_ticks)
            ax.set_yticklabels(y_ticks + 1, rotation=45)
            
            # Save the heatmap as a .png file with DPI 300
            output_file_path = os.path.join(output_dir, f"{filename.split('.')[0]}_heatmap.png")
            plt.savefig(output_file_path, dpi=300, bbox_inches='tight')
            
            # Close the current figure
            plt.close()
        except Exception as e:
            print(f"Error processing {filename}: {e}")
