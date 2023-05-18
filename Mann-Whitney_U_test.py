import os
import pandas as pd
from scipy.stats import mannwhitneyu

#--- Set Path
os.chdir(r"C:\Users\saySa\OneDrive\Desktop\test\betweenness")

file_names = ['2j4ddataresiduemap.xlsx', 
              '3umvdataresiduemap.xlsx', 
              '3zxsdataresiduemap.xlsx', 
              '4u63dataresiduemap.xlsx', 
              '4gu5dataresiduemap.xlsx']


for i in range(len(file_names)-1):
    for j in range(i+1, len(file_names)):
        group1_data = pd.read_excel(file_names[i])
        group2_data = pd.read_excel(file_names[j])

        #--- Extract the betweenness centrality values for each group
        group1_values = group1_data['Betweenness'].values
        group2_values = group2_data['Betweenness'].values

        #--- Perform the Mann-Whitney U test
        statistic, p_value = mannwhitneyu(group1_values, group2_values)

        print(f"Mann-Whitney U Test between {file_names[i]} and {file_names[j]}:")
        print(f"Test Statistic: {statistic:.4f}")
        print(f"P-value: {p_value:.4f}")

        #--- Check the significance based on the p-value
        alpha = 0.05  # Significance level
        if p_value > alpha:
            print("There is no significant difference between the two groups.")
        else:
            print("There is a significant difference between the two groups.")
        print()