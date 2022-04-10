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

#-- Set Path 
os.chdir("C:\\Users\\saySa\\OneDrive\\Desktop\\Multilayer_Network_Modelling\\Datasets\PCN_test\\dataset_1")

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

def fetch(file): 
    
    atom_type = MMCIF2Dict(file) ['_atom_site.group_PDB']
    atom_site = MMCIF2Dict(file) ['_atom_site.id']
    atom_symbol = MMCIF2Dict(file) ['_atom_site.type_symbol']
    atom_fullsymbol = MMCIF2Dict(file) ['_atom_site.label_atom_id']
    atom_altloc = MMCIF2Dict(file) ['_atom_site.label_alt_id']
    residue_name = MMCIF2Dict(file) ['_atom_site.label_comp_id']
    atom_assym = MMCIF2Dict(file) ['_atom_site.label_asym_id']
    residue_id = MMCIF2Dict(file) ['_atom_site.label_seq_id']
    insertion_code = MMCIF2Dict(file) ['_atom_site.pdbx_PDB_ins_code']
    x_coord = MMCIF2Dict(file) ['_atom_site.Cartn_x']
    y_coord = MMCIF2Dict(file) ['_atom_site.Cartn_y']
    z_coord = MMCIF2Dict(file) ['_atom_site.Cartn_z']
    occupancy = MMCIF2Dict(file) ['_atom_site.occupancy']
    b_factor = MMCIF2Dict(file) ['_atom_site.B_iso_or_equiv']
    model_number = MMCIF2Dict(file) ['_atom_site.pdbx_PDB_model_num']
    
    data_cif = pd.DataFrame(list(zip(atom_type, 
                    atom_site, 
                    atom_symbol, 
                    atom_fullsymbol, 
                    atom_altloc, 
                    residue_name, 
                    atom_assym,
                    residue_id, 
                    insertion_code,
                    x_coord, 
                    y_coord, 
                    z_coord, 
                    occupancy, 
                    b_factor, 
                    model_number)))
    
    #--- remove HETATM, TER
    data_cif = (data_cif[~data_cif[0].isin(['HETATM','TER'])])
    
    #---    remove extra models
    data_cif[14] = data_cif[14].astype(int)
    data_cif = data_cif.drop(data_cif[data_cif[14] >1].index)  
    
    #---   remove altloc 
    data_cif = data_cif [~data_cif [4].isin(['B','C','D','E'])]
    
    #---   remove insertion codes
    data_cif = data_cif [~data_cif [8].isin(['B','C','D','E'])]
    
    #--- Pick CA
    pick_CA =  data_cif[data_cif[3] == 'CA']
    
    #--- Pick XYZ coordinates
    xyz_coord = pick_CA[pick_CA.columns[9:12]]
    
    #--- Pick the residue numbers
    res =  pick_CA[7]
    
    #--- Pick atomsite IDs  
    atomid = pick_CA[1]
    
    #--- AssymID
    assymid = pick_CA[6]
    
    temp = pd.concat([res, xyz_coord, atomid, assymid], axis=1)
    return temp

def _network(file):
    temp = fetch(file)
    res = temp[temp.columns[4]]
    xyz_matrix = temp[temp.columns[1:4]]
    points = np.array(xyz_matrix)
    dist_condensed = pdist(points)

    prefix = "ATM"
    labels = [prefix + item for item in res]
    
    tuples =[]
    for item in combinations(labels,2):
        tuples.append(item)
        
    source = pd.DataFrame(tuples)[0]
    target = pd.DataFrame(tuples)[1]
        
    edgelist = pd.DataFrame({'source': source,
                              'target': target,
                              'weights': dist_condensed
                            })
    
    #--- Short Range Interaction Network (SRIN)
    #set cutoff for LRIN at 12

    
    label_index = list(range(1,len(labels)+1))   
    # label_index_dataframe = pd.DataFrame(zip(label_index, labels))
    k1 = pd.DataFrame(list(combinations(label_index,2)))
    k2 = pd.DataFrame(list(combinations(labels,2)))
    k2[2] = edgelist['weights']
    k2[3] = k1[0]
    k2[4] = k1[1]
    diff = k2[4] - k2[3]
    k2[5] = diff
     
    newlist1 = k2[k2[5]<=12]
    newlist2 = newlist1[newlist1[2]<=8]
    cutoff_list = pd.DataFrame(list(zip(newlist2[0],newlist2[1],newlist2[2])), 
                               columns = ['source', 'target','weights'])
    
 
     
    return cutoff_list


fetch('5jry.cif') 
temp1 = _network('5jry.cif')
        
    
def percolation(file):
    
    df_subset = {}    
    for i in range(1000):
        df_len = len(temp1)      
        remove_frac = int((i/1000)*df_len)
        indices = temp1.index.values
        drop_indices = np.random.choice((indices), remove_frac, replace=False)
        df_subset[i]  = temp1.drop(drop_indices)  
    return df_subset            
temp5 = percolation('5jry.cif')

def plot_percolation(file):
    
    ccs = []
    H = nx.Graph()
    for i in temp5:       
        H= nx.from_pandas_edgelist(temp5[i])
        nx.draw_shell(H)
        plt.show() 
        ccs.append(nx.connected_components(H))
     
    return ccs

fg = plot_percolation('5jry.cif')
dump = []
for i in fg:
    dump.append(list(i))
    
# count
gc_count = []
for i in range(len(dump)):
    gc_count.append(len((max(dump[i]))))
            
# plot
xdata = np.arange(0, 1, 0.001).tolist()
ydata= gc_count
plt.plot(xdata,ydata, color='black', linestyle='dashed', linewidth = 0.25,
          marker='o', markersize=2)
plt.title('5jry')
plt.xlabel("fraction of removed links")
plt.ylabel("size of the largest component")
