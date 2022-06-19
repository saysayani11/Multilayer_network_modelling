import os
import glob
import pandas as pd
import matplotlib.pyplot as plt

os.chdir(r"C:\Users\saySa\OneDrive\Desktop")

def read_textfile():
    with open('pdb_files_1.txt','r') as text:
        textfiles = text.readlines()
        suffix_cif = ".cif"
        pdbids=[]
        for element in textfiles:
            pdbids.append(element.strip())
        pdbids = ["{}{}".format(i,suffix_cif) for i in pdbids]
    return pdbids
cif_files = list(read_textfile())


for i in cif_files:
    
    var = pd.read_excel(i+".xlsx")	
    x = list(var['frac_of_removed_links'])
    y = list(var['size_of_largest_component_normalized'])
    plt.figure(figsize = (10,10), dpi = 600)
    
    
    plt.plot(x,y, color='blue', linestyle='dashed', linewidth = 0.2,
         marker='o', markersize=5)
    plt.title(i, size = 30)
    
    plt.xlabel("fraction of removed links", size = 20)
    plt.xticks(fontsize=16)
    plt.ylabel("size of the largest connected component (normalized)", size = 20)
    plt.yticks(fontsize=16)
    plt.show()