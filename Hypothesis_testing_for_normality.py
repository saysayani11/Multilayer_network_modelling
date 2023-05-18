import os
import pandas as pd
from scipy.stats import shapiro

#--- Set Path
os.chdir(r"C:\Users\saySa\OneDrive\Desktop\test\betweenness")

file_names = ['2j4ddataresiduemap.xlsx', 
              '3umvdataresiduemap.xlsx', 
              '3zxsdataresiduemap.xlsx', 
              '4u63dataresiduemap.xlsx', 
              '4gu5dataresiduemap.xlsx']

#--- Shapiro-Wilk test for each dataset
for file_name in file_names:
    data = pd.read_excel(file_name)
    betweenness_centrality = data['Betweenness'].values
    stat, p_value = shapiro(betweenness_centrality)
    print(f"Dataset: {file_name}")
    print(f"Test Statistic: {stat:.4f}")
    print(f"P-value: {p_value:.4f}")

    #--- Interpretation based on the p-value
    alpha = 0.05  #--- Significance level
    if p_value > alpha:
        print("The data appears to be normally distributed.")
    else:
        print("The data does not appear to be normally distributed.")
    print()
