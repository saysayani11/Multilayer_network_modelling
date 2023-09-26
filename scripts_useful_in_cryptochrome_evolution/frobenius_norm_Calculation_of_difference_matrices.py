import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Specify the directory containing the .csv subtraction matrices
directory = r"C:\Users\saySa\OneDrive\Desktop\task_1\Distance_matrix_Analysis\matrix_substraction\data"

frobenius_norms = []
matrix_numbers = []

for idx, filename in enumerate(os.listdir(directory)):
    if filename.endswith(".csv"):  # matrices are saved as .csv files
        file_path = os.path.join(directory, filename)
        matrix = pd.read_csv(file_path, header=None).values  # load matrix as numpy array
        frobenius_norm = np.linalg.norm(matrix, 'fro')
        frobenius_norms.append(frobenius_norm)
        matrix_numbers.append(idx+1)

# Create a DataFrame with the Frobenius norms data
df = pd.DataFrame({'Matrix Number': matrix_numbers, 'Frobenius Norm': frobenius_norms})

# Save the DataFrame to an Excel file
output_file_path = "frobenius_norms.xlsx"
df.to_excel(output_file_path, index=False)

# Plot the calculated Frobenius norms
plt.figure(figsize=(20,12))
plt.bar(matrix_numbers, frobenius_norms, color='b')
plt.xlabel('Matrix Number')
plt.ylabel('Frobenius Norm')
plt.title('Frobenius Norms of Difference Matrices')
plt.xticks(matrix_numbers, rotation=45, fontsize=5)  # rotate tick labels by 45 degrees and reduce font size
plt.tight_layout()
plt.savefig('frobenius_norms_plot.png', dpi=300)
plt.show()
