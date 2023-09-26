import pandas as pd
import os

# Set the working directory
os.chdir(r"C:\Users\saySa\OneDrive\Desktop\task_1\PDBeFOLD")

# Define the file paths
original_dataset_file_path = 'PDBeFOLD_alignment_residues_only.xlsx'
number_dataset_file_path = 'PDBeFOLD_closeness.xlsx'

# Read the datasets from Excel files
original_dataset_df = pd.read_excel(original_dataset_file_path, engine='openpyxl')
number_dataset_df = pd.read_excel(number_dataset_file_path, engine='openpyxl')

# Ensure the column names are the same
original_dataset_df.columns = number_dataset_df.columns

# Create a result DataFrame with the same structure as original_dataset_df
result_df = original_dataset_df.copy()

# Initialize an iterator for each column in number_dataset_df
iterators = {col: iter(number_dataset_df[col]) for col in number_dataset_df.columns}

# Iterate over the original_dataset_df and fill the values in result_df
for col in original_dataset_df.columns:
    for idx in original_dataset_df.index:
        # If the corresponding cell in original_dataset_df is not a space or NaN,
        # replace the value in result_df with the next value from the corresponding column in number_dataset_df
        if pd.notna(original_dataset_df.at[idx, col]) and original_dataset_df.at[idx, col] != ' ':
            result_df.at[idx, col] = next(iterators[col])

# Define the output file path
output_file_path = 'result_dataset_closeness.xlsx'

# Save the result DataFrame to an Excel file
result_df.to_excel(output_file_path, index=False, engine='openpyxl')
