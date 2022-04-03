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
    
    cutoff_list = edgelist[edgelist['weights']<8]
     
    #--- set coordinates for network
    chains =  list(set(temp[temp.columns[5]] ))
    count_chain_res = Counter(temp[temp.columns[5]])

    keys = labels
    dict_xyz = (xyz_matrix.astype(float).values.tolist())
    dict_network = dict(list((zip(keys, dict_xyz))))


    #--- Plot the Network 
    H = nx.Graph()
    H = nx.from_pandas_edgelist(cutoff_list)
    pos = dict_network
    temp = list(H.edges())

    edge_xyz = [(pos[u], pos[v]) for u, v in H.edges()]
    t = np.asarray(edge_xyz)

    fig = plt.figure(figsize=(20,20))
    ax = fig.add_subplot(111, projection="3d")
        
    xs = np.array(xyz_matrix[9],dtype=float)
    ys = np.array(xyz_matrix[10],dtype=float)
    zs = np.array(xyz_matrix[11],dtype=float)

    ax.scatter(xs,ys,zs)
    
    for vizedge in t:
        ax.plot(*vizedge.T, linewidth=0.8,color="red")
    plt.title(str(file))  
    
    return cutoff_list


fetch('5jry.cif') 
temp1 = _network('5jry.cif')
        
def frac(n):
    df_len = len(temp1)
    remove_frac_1 = int(n*df_len)
    indices = temp1.index.values
    drop_indices = np.random.choice((indices), remove_frac_1, replace=False)
    df_subset = temp1.drop(drop_indices) 
    return df_subset
    
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


















#==============================================================================

# for i in cif_files:  
#     fetch(i)
#     _network(i)
#     percolation(i)
#     plot_percolation(i)

#==============================================================================
   

    # remove_frac_1 = int(0.1*df_len)
    # indices = temp1.index.values
    # drop_indices = np.random.choice((indices), remove_frac_1, replace=False)
    # df_subset = temp1.drop(drop_indices)   
    
    # H = nx.Graph()
    # H = nx.from_pandas_edgelist(df_subset)

    # fraclist = np.arange(0, 1, 0.01).tolist()
    # remove = np.multiply(fraclist,df_len)
    # indices = temp1.index.values

    # fraclist = np.arange(0, 1.01, 0.01).tolist()
    # remove = np.floor(np.multiply(fraclist,df_len))[::-1]
    
    # select = {}
   
    # -- remove k nodes  

 