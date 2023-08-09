import os
import numpy as np
import pandas as pd
import networkx as nx
from scipy.spatial.distance import pdist
import matplotlib.pyplot as plt
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from itertools import combinations
import seaborn as sns

#-- Set Path
os.chdir(r"C:\Users\saySa\OneDrive\Desktop\task_1\PSNs\networks")
plt.rcParams.update({'figure.max_open_warning': 0, 'figure.dpi': 300})

def read_textfile():
    with open('structures_files.txt','r') as text:
        textfiles = text.readlines()
        suffix_cif = ".cif"
        pdbids=[]
        for element in textfiles:
            pdbids.append(element.strip())
        pdbids = ["{}{}".format(i,suffix_cif) for i in pdbids]
    return pdbids

def fetch(file):
    mmcif_dict = MMCIF2Dict(file)
    data_cif = pd.DataFrame(list(zip(
        mmcif_dict['_atom_site.group_PDB'],
        mmcif_dict['_atom_site.id'],
        mmcif_dict['_atom_site.type_symbol'],
        mmcif_dict['_atom_site.label_atom_id'],
        mmcif_dict['_atom_site.label_alt_id'],
        mmcif_dict['_atom_site.label_comp_id'],
        mmcif_dict['_atom_site.label_asym_id'],
        mmcif_dict['_atom_site.label_seq_id'],
        mmcif_dict['_atom_site.pdbx_PDB_ins_code'],
        mmcif_dict['_atom_site.Cartn_x'],
        mmcif_dict['_atom_site.Cartn_y'],
        mmcif_dict['_atom_site.Cartn_z'],
        mmcif_dict['_atom_site.occupancy'],
        mmcif_dict['_atom_site.B_iso_or_equiv'],
        mmcif_dict['_atom_site.pdbx_PDB_model_num']
    )))
   
    data_cif = (data_cif[~data_cif[0].isin(['HETATM','TER'])])
    data_cif[14] = data_cif[14].astype(int)
    data_cif = data_cif.drop(data_cif[data_cif[14] >1].index)
    data_cif = data_cif[~data_cif[4].isin(['B','C','D','E'])]
    data_cif = data_cif[~data_cif[8].isin(['B','C','D','E'])]
   
    pick_CA =  data_cif[data_cif[3] == 'CA']
    xyz_coord = pick_CA[pick_CA.columns[9:12]]
    res = pick_CA[7]
    atomid = pick_CA[1]
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
    
    tuples = []
    for item in combinations(labels, 2):
        tuples.append(item)
        
    edgelist = pd.DataFrame({
        'source': [t[0] for t in tuples],
        'target': [t[1] for t in tuples],
        'distance': dist_condensed
    })
    
    cutoff_list = edgelist[edgelist['distance']<8]
    return cutoff_list



def _degdist2(cif_files):
    # Initialize a list to store degree frequencies for each file
    degree_freqs = []

    for file in cif_files:
        temp = _network(file)
        H = nx.Graph()
        H = nx.from_pandas_edgelist(temp)
        degree_freq = nx.degree_histogram(H)
        degree_freqs.append(degree_freq)
    
    n_files = len(cif_files)
    color_palette = sns.color_palette("flare", n_files)

    plt.figure(figsize=(14, 8), dpi=300)
    
    for i, degree_freq in enumerate(degree_freqs):
        label = os.path.splitext(cif_files[i])[0]
        if label == "bcpd2_cut":
            label = "bCPD2"
        x_values = list(range(len(degree_freq)))
        plt.scatter(x_values, degree_freq, label=label, color=color_palette[i])

    plt.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=14)
    plt.title("Degree Distribution for all 18 Structures", fontsize=18)
    plt.xlabel("Degree", fontsize=16)
    plt.ylabel("Frequency", fontsize=16)
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.tight_layout()
    plt.show()

cif_files = list(read_textfile())
_degdist2(cif_files)

