import os
import random
import numpy as np
import pandas as pd
import networkx as nx
from networkx.generators.degree_seq import expected_degree_graph

import matplotlib.pyplot as plt
from matplotlib import cm
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from scipy.spatial.distance import pdist, squareform
from itertools import combinations
from collections import Counter



plt.rcParams.update({'figure.max_open_warning': 0})

#-- Set Path
os.chdir("/mnt/808e843c-fb08-404a-9f45-efd817d81a77/Sayantoni/networks/cif_files_1758/")

def read_textfile():
    with open('1758_dataset_19_04_2022','r') as text:
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
    points = np.array(xyz_matrix).astype(float)
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
                              'distance': dist_condensed
                            })
    
    filename = str(file)
    cutoff_list = edgelist[edgelist['distance']<8]
    
    
    return cutoff_list

def _degdist2(file):
    temp = _network(file)
    H = nx.Graph()
    H = nx.from_pandas_edgelist(temp)
    degree_freq = nx.degree_histogram(H)
    degrees = range(len(degree_freq))
    
    plt.figure(figsize=(12, 8)) 
    plt.loglog(degrees, degree_freq,'go-', label = "original network") 

    
    print ("degree seqence for original network = " , list(degree_freq))

    
    ############################ERDOS RENYI NETWORK WITH NO DEGREE PRESERVATION (full randomization)
    n = nx.number_of_nodes(H)
    m = nx.number_of_edges(H)
    K = nx.Graph()
    K = nx.gnm_random_graph(n,m)
    
    degree_freq_er = nx.degree_histogram(K)
    degrees_er = range(len(degree_freq_er))

    print ("degree seqence for er1 = " , list(degree_freq_er))
    plt.loglog(degrees_er, degree_freq_er,'bo-',label = 'fully randomized') 

    # #############################ERDOS RENYI NETWORK WITH EXPECTED DEGREE 
    # G2 = expected_degree_graph([deg for (k, deg) in H.degree()], selfloops=False)
    # degree_freq_er2 = nx.degree_histogram(G2)
    # degrees_er2 = range(len(degree_freq_er2))
    
    # print ("degree seqence for er2 = " , list(degree_freq_er2))
    # plt.loglog(degrees_er2, degree_freq_er2,'ro-') 
    
    
    #############################ERDOS RENYI NETWORK WITH DEGREE PRESERVING RANDOMIZATION
    degrees = [val for (node, val) in H.degree()]
    k = random.sample(degrees, len(degrees))
    G3 = nx.configuration_model(k)
    G3.remove_edges_from(nx.selfloop_edges(G3))
    degree_freq_n = nx.degree_histogram(G3)
    degrees_n = [val for (node, val) in G3.degree()]
    degree_counts_n = Counter(degrees_n)   
    x1, y1 = zip(*degree_counts_n.items())   
    dictionary1 = dict(zip(x1, y1))
    sorted_dict1 = dict(sorted(dictionary1.items()))
    
    d1 = list(sorted_dict1.keys())
    f1 = list(sorted_dict1.values())
    plt.plot(d1,f1, 'k+--',label = 'degree preserving') 
    
    print ("degree seqence for randomize network = " , list(degree_freq_n))

    plt.xlabel("node degree")
    plt.ylabel("frequency")
    # plt.show()
    
    #--- SAVE PLOT
    filename = str(file) [:-4]
    directory = "/mnt/808e843c-fb08-404a-9f45-efd817d81a77/Sayantoni/networks/results/" + filename
    subdirectory = directory + "/degree_distribution"
    os.chdir(subdirectory) 
    
    plt.title("Degree Distribution")

    leg = plt.legend()
    

    plt.legend(loc='upper left', frameon=False)
    plt.savefig(filename + ".tiff", dpi = 300)
    

    
    #--- SAVE DATA
    df = pd.DataFrame (zip(degrees, degree_freq, degree_freq_er, degrees_er, d1,f1))
    df.to_excel(filename + "degree_dist.xlsx")
    os.chdir("/mnt/808e843c-fb08-404a-9f45-efd817d81a77/Sayantoni/networks/cif_files_1758")
    

for i in cif_files:
    _network(i)
    # _degdist1(i)
    _degdist2(i)
