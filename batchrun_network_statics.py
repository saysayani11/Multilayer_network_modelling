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

def _distmatrix(file):
    temp = fetch(file)
    xyz_matrix = temp[temp.columns[1:4]]
    points = np.array(xyz_matrix)
    dist_condensed = pdist(points)
    dist = squareform(dist_condensed)
    return dist

def _adjacencymatrix(file):
    adj_matrix = _distmatrix(file)>=8
    np.fill_diagonal(adj_matrix,1)
    #--- Plot Adjacency Matrix
    plt.matshow(adj_matrix,cmap=cm.binary)
    plt.title(str(file))   
    return adj_matrix

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

ty = _network('2pdd.cif')

def _degdist(file):
    temp = _network(file)
    H = nx.Graph()
    H = nx.from_pandas_edgelist(temp)
    #-- degrees of each node
    degrees = [val for (node, val) in H.degree()]
    degree_counts = Counter(degrees)                                                                                                 
    x, y = zip(*degree_counts.items())   
    
    plt.figure(figsize = (12,8))
    plt.xlabel('degree')                                                                                                             
                                                                                                               
    plt.xlim(0, max(x)+10)  
                                                                                                           
    plt.ylabel('frequency')                                                                                                          
                                                                                                            
    plt.ylim(0, max(y)+10)                                                                                                             
                                                                                                                                     # do plot                                                                                                                        
    plt.scatter(x, y, marker='o')                                                                                                    
    plt.show()
    
    return x,y

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



# def _BCA(file):
#     temp = _network(file)
#     H = nx.Graph()
#     H = nx.from_pandas_edgelist(temp)
#     bca = nx.betweenness_centrality(H)  
#     dpos = list(bca.values())
    
#     temp = fetch(file)
#     res = temp[temp.columns[4]]
#     xyz_matrix = temp[temp.columns[1:4]]
    
#     prefix = "ATM"
#     labels = [prefix + item for item in res]
#     keys = labels
#     dict_xyz = (xyz_matrix.astype(float).values.tolist())
#     dict_network = dict(list((zip(keys, dict_xyz))))
    
#     fig = plt.figure(figsize=(20,20))   
#     ax = fig.add_subplot(111, projection="3d")
        
#     xs = np.array(xyz_matrix[9],dtype=float)
#     ys = np.array(xyz_matrix[10],dtype=float)
#     zs = np.array(xyz_matrix[11],dtype=float) 
#     t = ax.scatter(xs,ys,zs, s = 90, c = dpos, cmap = 'gist_yarg')
#     fig.colorbar(t, ax=ax)
#     plt.title(str(file))    
        
#     pos = dict_network
#     temp = list(H.edges())

#     edge_xyz = [(pos[u], pos[v]) for u, v in H.edges()]
#     t = np.asarray(edge_xyz)
#     for vizedge in t:
#         ax.plot(*vizedge.T, linewidth=0.4,color="black")
#     # print("Betweenness centrality of " + str(file) + ": ", bca)
#     return bca

# def _CC(file):
#     temp = _network(file)
#     H = nx.Graph()
#     H = nx.from_pandas_edgelist(temp)
#     cc = nx.closeness_centrality(H)  
#     dpos = list(cc.values())
    
#     temp = fetch(file)
#     res = temp[temp.columns[4]]
#     xyz_matrix = temp[temp.columns[1:4]]
    
#     prefix = "ATM"
#     labels = [prefix + item for item in res]
#     keys = labels
#     dict_xyz = (xyz_matrix.astype(float).values.tolist())
#     dict_network = dict(list((zip(keys, dict_xyz))))
    
#     fig = plt.figure(figsize=(20,20))   
#     ax = fig.add_subplot(111, projection="3d")
        
#     xs = np.array(xyz_matrix[9],dtype=float)
#     ys = np.array(xyz_matrix[10],dtype=float)
#     zs = np.array(xyz_matrix[11],dtype=float)
#     t = ax.scatter(xs,ys,zs, s = 100, c = dpos, cmap = cm.gist_yarg)
#     fig.colorbar(t, ax=ax)
#     plt.title(str(file))    
        
#     pos = dict_network
#     temp = list(H.edges())

#     edge_xyz = [(pos[u], pos[v]) for u, v in H.edges()]
#     t = np.asarray(edge_xyz)
#     for vizedge in t:
#         ax.plot(*vizedge.T, linewidth=0.4,color="black")

#     # print("Closeness Centrality of " + str(file) + ": ", cc)
#     return cc
    
# def _ECA(file):
#     temp = _network(file)
#     H = nx.Graph()
#     H = nx.from_pandas_edgelist(temp)
#     eca = nx.eigenvector_centrality_numpy(H)  
#     dpos = list(eca.values())
    
#     temp = fetch(file)
#     res = temp[temp.columns[4]]
#     xyz_matrix = temp[temp.columns[1:4]]
    
#     prefix = "ATM"
#     labels = [prefix + item for item in res]
#     keys = labels
#     dict_xyz = (xyz_matrix.astype(float).values.tolist())
#     dict_network = dict(list((zip(keys, dict_xyz))))
    
#     fig = plt.figure(figsize=(20,20))   
#     ax = fig.add_subplot(111, projection="3d")
        
#     xs = np.array(xyz_matrix[9],dtype=float)
#     ys = np.array(xyz_matrix[10],dtype=float)
#     zs = np.array(xyz_matrix[11],dtype=float)
#     t = ax.scatter(xs,ys,zs, s = 100, c = dpos, cmap = cm.gist_yarg)
#     fig.colorbar(t, ax=ax)
#     plt.title(str(file))    
        
#     pos = dict_network
#     temp = list(H.edges())

#     edge_xyz = [(pos[u], pos[v]) for u, v in H.edges()]
#     t = np.asarray(edge_xyz)
#     for vizedge in t:
#         ax.plot(*vizedge.T, linewidth=0.4,color="black")

#     # print("eigenvector centrality of " + str(file) + ": ", cc)
#     return eca

# def _CLA(file):
#     temp = _network(file)
#     H = nx.Graph()
#     H = nx.from_pandas_edgelist(temp)
#     cla = nx.clustering(H)  
#     dpos = list(cla.values())
    
#     temp = fetch(file)
#     res = temp[temp.columns[4]]
#     xyz_matrix = temp[temp.columns[1:4]]
    
#     prefix = "ATM"
#     labels = [prefix + item for item in res]
#     keys = labels
#     dict_xyz = (xyz_matrix.astype(float).values.tolist())
#     dict_network = dict(list((zip(keys, dict_xyz))))
    
#     fig = plt.figure(figsize=(20,20))   
#     ax = fig.add_subplot(111, projection="3d")
        
#     xs = np.array(xyz_matrix[9],dtype=float)
#     ys = np.array(xyz_matrix[10],dtype=float)
#     zs = np.array(xyz_matrix[11],dtype=float)
#     t = ax.scatter(xs,ys,zs, s = 100, c = dpos, cmap = cm.gist_yarg)
#     fig.colorbar(t, ax=ax)
#     plt.title(str(file))    
        
#     pos = dict_network
#     temp = list(H.edges())

#     edge_xyz = [(pos[u], pos[v]) for u, v in H.edges()]
#     t = np.asarray(edge_xyz)
#     for vizedge in t:
#         ax.plot(*vizedge.T, linewidth=0.4,color="black")

#     # print("eigenvector centrality of " + str(file) + ": ", cc)
#     return cla
    
# def _is_small_world(file):  
#     H = nx.Graph()
#     H = nx.from_pandas_edgelist(_network(file))
#     sigma_val = nx.sigma(H, niter=10, nrand=10, seed=None)
#     return sigma_val
    
# # ---- WRITE DATA ----
# def _write(file):
#     df = pd.DataFrame.from_dict(_DCA(file), orient = 'index')
#     with pd.ExcelWriter('Dataset1.xlsx', engine='openpyxl', mode='a') as writer: 
#           df.to_excel(writer)     
#     # df.to_excel('Dataset1.xlsx', mode = 'a')
#     return df



fetch('2pdd.cif')
_distmatrix('2pdd.cif')
_adjacencymatrix('2pdd.cif')
_degdist('2pdd.cif')
    # _DCA(i)
    # _BCA(i)
    # _CC(i)
    # _ECA(i)
    # _CLA(i)
    # _is_small_world(i)
    
    
    
    
    



  






