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

    
def percolation(file):
    
    df_subset1 = {}   
    df_subset2 = {}
    df_subset3 = {}
    df_subset4 = {}
    df_subset5 = {}
    df_subset6 = {}   
    df_subset7 = {}
    df_subset8 = {}
    df_subset9 = {}
    df_subset10 = {}
    
    
    for i in range(1000):
        df_len = len(temp1)      
        remove_frac = int((i/1000)*df_len)
        indices = temp1.index.values
        drop_indices1 = np.random.choice((indices), remove_frac, replace=False)
        drop_indices2 = np.random.choice((indices), remove_frac, replace=False)
        drop_indices3 = np.random.choice((indices), remove_frac, replace=False)
        drop_indices4 = np.random.choice((indices), remove_frac, replace=False)
        drop_indices5 = np.random.choice((indices), remove_frac, replace=False)
        drop_indices6 = np.random.choice((indices), remove_frac, replace=False)
        drop_indices7 = np.random.choice((indices), remove_frac, replace=False)
        drop_indices8 = np.random.choice((indices), remove_frac, replace=False)
        drop_indices9 = np.random.choice((indices), remove_frac, replace=False)
        drop_indices10 = np.random.choice((indices), remove_frac, replace=False)
        
        df_subset1[i]  = temp1.drop(drop_indices1)  
        df_subset2[i]  = temp1.drop(drop_indices2) 
        df_subset3[i]  = temp1.drop(drop_indices3) 
        df_subset4[i]  = temp1.drop(drop_indices4) 
        df_subset5[i]  = temp1.drop(drop_indices5) 
        df_subset6[i]  = temp1.drop(drop_indices6)  
        df_subset7[i]  = temp1.drop(drop_indices7) 
        df_subset8[i]  = temp1.drop(drop_indices8) 
        df_subset9[i]  = temp1.drop(drop_indices9) 
        df_subset10[i]  = temp1.drop(drop_indices10) 
        
    return df_subset1, df_subset2, df_subset3, df_subset4, df_subset5, df_subset6, df_subset7, df_subset8, df_subset9, df_subset10        
temp5 = percolation('5jry.cif')

def plot_percolation(file):
    
    ccs1 = [] 
    ccs2 = []
    ccs3 = []
    ccs4 = []
    ccs5 = []
    ccs6 = [] 
    ccs7 = []
    ccs8 = []
    ccs9 = []
    ccs10 = []
    
    H1 = nx.Graph()
    H2 = nx.Graph()
    H3 = nx.Graph()
    H4 = nx.Graph()
    H5 = nx.Graph()
    H6 = nx.Graph()
    H7 = nx.Graph()
    H8 = nx.Graph()
    H9 = nx.Graph()
    H10 = nx.Graph()

    df_temp5_0 = temp5[0]
    for i in df_temp5_0:       
        H1 = nx.from_pandas_edgelist(df_temp5_0[i])
        ccs1.append(nx.connected_components(H1))
        
    df_temp5_1 = temp5[1]
    for i in df_temp5_1:       
        H2 = nx.from_pandas_edgelist(df_temp5_1[i])
        ccs2.append(nx.connected_components(H2))
    
    df_temp5_2 = temp5[2]
    for i in df_temp5_2:       
        H3 = nx.from_pandas_edgelist(df_temp5_2[i]) 
        ccs3.append(nx.connected_components(H3))
        
    df_temp5_3 = temp5[3]
    for i in df_temp5_3:       
        H4 = nx.from_pandas_edgelist(df_temp5_3[i])
        ccs4.append(nx.connected_components(H4))
                           
    df_temp5_4 = temp5[4]
    for i in df_temp5_4:       
        H5 = nx.from_pandas_edgelist(df_temp5_4[i])
        ccs5.append(nx.connected_components(H5))
        
    df_temp5_5 = temp5[5]
    for i in df_temp5_5:       
        H6 = nx.from_pandas_edgelist(df_temp5_5[i])
        ccs6.append(nx.connected_components(H6))
        
    df_temp5_6 = temp5[6]
    for i in df_temp5_6:       
        H7 = nx.from_pandas_edgelist(df_temp5_6[i])
        ccs7.append(nx.connected_components(H7))
    
    df_temp5_7 = temp5[7]
    for i in df_temp5_7:       
        H8 = nx.from_pandas_edgelist(df_temp5_7[i]) 
        ccs8.append(nx.connected_components(H8))
        
    df_temp5_8 = temp5[8]
    for i in df_temp5_8:       
        H9 = nx.from_pandas_edgelist(df_temp5_8[i])
        ccs9.append(nx.connected_components(H9))
                           
    df_temp5_9 = temp5[9]
    for i in df_temp5_9:       
        H10 = nx.from_pandas_edgelist(df_temp5_9[i])
        ccs10.append(nx.connected_components(H10))
     
    return ccs1,ccs2,ccs3,ccs4,ccs5,ccs6,ccs7,ccs8,ccs9,ccs10


fg = plot_percolation('5jry.cif')
fg1 = fg[0]
fg2 = fg[1]
fg3 = fg[2]
fg4 = fg[3]
fg5 = fg[4]
fg6 = fg[5]
fg7 = fg[6]
fg8 = fg[7]
fg9 = fg[8]
fg10 = fg[9]

dump1 = []
for i in fg1:
    dump1.append(list(i))
    
dump2 = []
for i in fg2:
    dump2.append(list(i))
    
dump3 = []
for i in fg3:
    dump3.append(list(i))
    
dump4 = []
for i in fg4:
    dump4.append(list(i))

dump5 = []
for i in fg5:
    dump5.append(list(i))
    
dump6 = []
for i in fg6:
    dump6.append(list(i))
    
dump7 = []
for i in fg7:
    dump7.append(list(i))
    
dump8 = []
for i in fg8:
    dump8.append(list(i))
    
dump9 = []
for i in fg9:
    dump9.append(list(i))

dump10 = []
for i in fg10:
    dump10.append(list(i))
    

    
#--- count
gc_count1 = []
for i in range(len(dump1)):
    gc_count1.append(len((max(dump1[i]))))
    
gc_count2 = []
for i in range(len(dump2)):
    gc_count2.append(len((max(dump2[i]))))
    
gc_count3 = []
for i in range(len(dump3)):
    gc_count3.append(len((max(dump3[i]))))
    
gc_count4 = []
for i in range(len(dump4)):
    gc_count4.append(len((max(dump4[i]))))

gc_count5 = []
for i in range(len(dump5)):
    gc_count5.append(len((max(dump5[i]))))
    
gc_count6 = []
for i in range(len(dump6)):
    gc_count6.append(len((max(dump6[i]))))
    
gc_count7 = []
for i in range(len(dump7)):
    gc_count7.append(len((max(dump7[i]))))
    
gc_count8 = []
for i in range(len(dump8)):
    gc_count8.append(len((max(dump8[i]))))
    
gc_count9 = []
for i in range(len(dump9)):
    gc_count9.append(len((max(dump9[i]))))

gc_count10 = []
for i in range(len(dump10)):
    gc_count10.append(len((max(dump10[i]))))
    
    
#--- plot
xdata = np.arange(0, 1, 0.001).tolist()
ydata1 = gc_count1
ydata2 = gc_count2
ydata3 = gc_count3
ydata4 = gc_count4
ydata5 = gc_count5
ydata6 = gc_count6
ydata7 = gc_count7
ydata8 = gc_count8
ydata9 = gc_count9
ydata10 = gc_count10

data_tuples = list(zip(xdata,ydata1,ydata2,ydata3,ydata4,ydata5,ydata6,ydata7,
                                     ydata8,ydata9,ydata10))

percolation_dataframe = pd.DataFrame(data_tuples)


#--- PLOT 

plt.scatter(xdata, ydata1, color='black', s = 2)
plt.scatter(xdata, ydata2, color='black', s = 2)
plt.scatter(xdata, ydata3, color='black', s = 2)
plt.plot(xdata, ydata2, color='black', markersize=2)
plt.plot(xdata, ydata3, color='black', markersize=2)
plt.plot(xdata, ydata4, color='black', markersize=2)
plt.plot(xdata, ydata5, color='black', markersize=2)
plt.xlabel("fraction of removed links")
plt.ylabel("size of the largest component")
plt.show()

percolation_dataframe.to_excel('perc_replicate_1.xlsx')




         
