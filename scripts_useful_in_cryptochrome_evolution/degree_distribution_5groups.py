from __future__ import print_function
import os
import random
import numpy as np
import pandas as pd
import networkx as nx
from networkx.generators.degree_seq import expected_degree_graph
from sklearn.preprocessing import normalize
from sklearn import preprocessing
import matplotlib.pyplot as plt
from matplotlib import cm
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from scipy.spatial.distance import pdist, squareform
from itertools import combinations
from collections import Counter

plt.rcParams.update({'figure.max_open_warning': 0, 'figure.dpi': 300})

#-- Set Path
os.chdir(r"C:\Users\saySa\OneDrive\Desktop\task_1\PSNs\networks")

def read_textfile():
    with open('structures_files.txt','r') as text:
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

import seaborn as sns

# Define your groups
group1 = ['3ZXS.cif', '4DJA.cif']
group2 = ['3UMV.cif', '2XRY.cif', 'bcpd2_cut.cif']
group3 = ['4GU5.cif', '5ZM0.cif', '7AYV.cif', '6PTZ.cif', '5T5X.cif', '7V8Y.cif']
group4 = ['1NP7.cif', '2J4D.cif']
group5 = ['1IQR.cif', '1DNP.cif', '4U63.cif', '1U3C.cif', '6K8I.cif']

# Define your groups list
groups = [group1, group2, group3, group4, group5]
group_names = ["Group 1 (6-4 like /FeS)", "Group 2 (CPD 2)", "Group 3 (3.6-4/Animal)", "Group 4 (DASH)", "Group 5 (5.CPD1/3)"]

# Create a subplot for each group
fig, axs = plt.subplots(3, 2, figsize=(18, 24))  # Adjusted size for better display

# Create a colormap
n_files = len(max(groups, key=len))  # Number of colors will be equal to the size of the largest group
cm = sns.color_palette("rocket_r", n_files, as_cmap=True)


# Create a subplot for each group
fig, axs = plt.subplots(1, len(groups), figsize=(18 * len(groups), 12))
plt.rcParams.update({'font.size': 50})

# Generate degree distribution for each group
for idx, (group, group_name) in enumerate(zip(groups, group_names)):
    ax = axs[idx]
    degree_freqs = []

    for file in group:
        temp = _network(file)
        H = nx.Graph()
        H = nx.from_pandas_edgelist(temp)
        degree_freq = nx.degree_histogram(H)

        degree_freqs.append(degree_freq)
    
    for i, degree_freq in enumerate(degree_freqs):
        label = os.path.splitext(group[i])[0]
        if label == "bcpd2_cut":
            label = "bCPD2"
        ax.plot(degree_freq, 'o', label=label, color=cm(i/n_files), markersize=15)  # Modification: Plot only dots without lines

    ax.tick_params(axis='both', which='major', labelsize=30)
    ax.tick_params(axis='both', which='minor', labelsize=30)
    ax.legend(prop={'size': 50})

plt.tight_layout()
plt.show()
