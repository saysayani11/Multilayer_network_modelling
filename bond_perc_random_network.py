import os
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib import cm
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from scipy.spatial.distance import pdist, squareform
from itertools import combinations
from collections import Counter
plt.rcParams.update({'figure.max_open_warning': 0})

n = 477
m = 2573

import random
random_sourcelist = []
for i in range(m):
     random_sourcelist.append(random.randint(1,477))
     
random_targetlist = []
for i in range(m):
     random_targetlist.append(random.randint(1,477))
     
temp = list(zip(random_sourcelist, random_targetlist))    
edgelist = pd.DataFrame(temp, columns = ['source', 'target'])
    
df_subset = {}    
for i in range(1000):
    df_len = len(edgelist)      
    remove_frac = int((i/1000)*df_len)
    indices = edgelist.index.values
    drop_indices = np.random.choice((indices), remove_frac, replace=False)
    df_subset[i]  = edgelist.drop(drop_indices)  
    

ccs = []
H = nx.Graph()
for i in df_subset:       
    H= nx.from_pandas_edgelist(df_subset[i])
    nx.draw_shell(H)
    plt.show() 
    ccs.append(nx.connected_components(H))
    
dump = []
for i in ccs:
    dump.append(list(i))
    
gc_count = []
for i in range(len(dump)):
    gc_count.append(len((max(dump[i]))))
            
    
xdata = np.arange(0, 1, 0.001).tolist()
ydata= (gc_count - (np.min(gc_count))) / (np.max(gc_count) - np.min(gc_count))
plt.scatter(xdata,ydata, color='black')


plt.title('5jry')
plt.xlabel("fraction of removed links")
plt.ylabel("size of the largest component")
