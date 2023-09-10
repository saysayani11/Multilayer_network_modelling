import pandas as pd
import os

# Set working directory and read the Excel file
os.chdir(r"C:\Users\saySa\OneDrive\Desktop\task_1\MSA\consurf_data")
file_path = 'CONSURF_RANK.xlsx'
df = pd.read_excel(file_path)

# Replace asterisks and convert to float type
df = df.applymap(lambda x: x.replace('*', '') if isinstance(x, str) else x)
df = df.astype(float)

# Count values for each category
variable_residue_count = df.apply(lambda col: (col.between(1, 4)).sum())
average_conservation_residue_count = df.apply(lambda col: (col == 5).sum())
conserved_residue_count = df.apply(lambda col: (col.between(6, 9)).sum())

# Create a summary DataFrame
summary_df = pd.DataFrame({
    'Variable_Residue_Count': variable_residue_count,
    'Average_Conservation_Residue_Count': average_conservation_residue_count,
    'Conserved_Residue_Count': conserved_residue_count
})

# Save the summary to an Excel file
summary_df.to_excel('summary_counts.xlsx', index=True)
