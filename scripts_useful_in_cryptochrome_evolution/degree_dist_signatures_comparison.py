import os
import numpy as np
import pandas as pd
from scipy import stats
from scipy.integrate import simps
from scipy.stats import gaussian_kde
from Bio.PDB import MMCIFParser, Superimposer, PDBIO
from Bio.PDB.PDBIO import PDBIO
import networkx as nx
from Bio.PDB import PDBIO
import matplotlib.pyplot as plt
from matplotlib import cm
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.MMCIFParser import MMCIFParser
from scipy.spatial.distance import pdist, squareform
import itertools
from Bio.PDB.Superimposer import Superimposer
from itertools import combinations
from collections import Counter
import seaborn as sns
plt.rcParams.update({'figure.max_open_warning': 0, 'figure.dpi': 300})

os.chdir(r"C:\Users\saySa\OneDrive\Desktop\task_1\PSNs\networks")

def read_textfile():
    with open('test.txt','r') as text:
        textfiles = text.readlines()
        suffix_cif = ".cif"
        pdbids=[]
        for element in textfiles:
            pdbids.append(element.strip())
        pdbids = ["{}{}".format(i,suffix_cif) for i in pdbids]
    return pdbids
cif_files = list(read_textfile())

def fetch(file):
    mmcif_dict = MMCIF2Dict(file)
    data_cif = pd.DataFrame({
        'atom_type': mmcif_dict['_atom_site.group_PDB'],
        'atom_site': mmcif_dict['_atom_site.id'],
        'atom_symbol': mmcif_dict['_atom_site.type_symbol'],
        'atom_fullsymbol': mmcif_dict['_atom_site.label_atom_id'],
        'atom_altloc': mmcif_dict['_atom_site.label_alt_id'],
        'residue_name': mmcif_dict['_atom_site.label_comp_id'],
        'atom_assym': mmcif_dict['_atom_site.label_asym_id'],
        'residue_id': mmcif_dict['_atom_site.label_seq_id'],
        'insertion_code': mmcif_dict['_atom_site.pdbx_PDB_ins_code'],
        'x_coord': mmcif_dict['_atom_site.Cartn_x'],
        'y_coord': mmcif_dict['_atom_site.Cartn_y'],
        'z_coord': mmcif_dict['_atom_site.Cartn_z'],
        'occupancy': mmcif_dict['_atom_site.occupancy'],
        'b_factor': mmcif_dict['_atom_site.B_iso_or_equiv'],
        'model_number': mmcif_dict['_atom_site.pdbx_PDB_model_num']
    })
    
    #-- Removal / Filter
    data_cif = data_cif[
            ~data_cif['atom_type'].isin(['HETATM', 'TER']) &
            (data_cif['model_number'].astype(int) == 1) &
            ~data_cif['atom_altloc'].isin(['B', 'C', 'D', 'E']) &
            ~data_cif['insertion_code'].isin(['B', 'C', 'D', 'E']) &
            (data_cif['atom_assym'] == 'A') &
            (data_cif['atom_fullsymbol'] == 'CA')
                       ]
    
    data_cif = data_cif[['residue_id', 'x_coord', 'y_coord', 'z_coord', 'atom_site', 'atom_assym', 'residue_name']]
    data_cif.columns = ['residue_no', 'x_coord', 'y_coord', 'z_coord', 'ca_atom_number', 'chain', 'residue_name']
    
    return data_cif

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
    cutoff_list = edgelist[edgelist['distance']<8]
    return cutoff_list

def _degree_distribution(file):
    # Create a network using the existing _network function
    edge_list = _network(file)
    
    # Create a graph using NetworkX
    G = nx.from_pandas_edgelist(edge_list, 'source', 'target', ['distance'])
    
    # Calculate the degree for each node
    degrees = dict(G.degree())
    
    # Convert the dictionary to a pandas DataFrame
    degree_df = pd.DataFrame(list(degrees.items()), columns=['Node', 'Degree'])
    
    # Optionally, store the degree distribution in a CSV file
    # This might be helpful if you wish to save the degree distribution of each protein structure in separate files
    # degree_df.to_csv(file.split('.')[0] + '_degree_distribution.csv', index=False)
    
    return degree_df

# Define a function to plot degree distributions for a list of protein files
def plot_degree_distributions(protein_files, figure_title, save_name):
    plt.figure(figsize=(10, 6))
    for cif_file in protein_files:
        degree_df = _degree_distribution(cif_file)
        # Use seaborn to plot the KDE
        sns.kdeplot(degree_df['Degree'], label=cif_file.split('.')[0], shade=True)
    
    plt.title(figure_title)
    plt.xlabel('Degree')
    plt.ylabel('Density')
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_name, dpi=300)
    plt.close()

# Define sets of proteins for comparison
protein_sets = [
    (["7ayv.cif", "5zm0.cif"], "Degree Distribution: 7ayv vs 5zm0", "7ayv_5zm0.png"),
    (["3umv.cif", "2xry.cif", "bcpd2_cut.cif"], "Degree Distribution: 3umv, 2xry vs bacterial CPD2", "3umv_2xry_bacterial_cpd2.png"),
    (["3zxs.cif", "4dja.cif"], "Degree Distribution: 3zxs vs 4dja", "3zxs_4dja.png"),
    (["2j4d.cif", "1np7.cif"], "Degree Distribution: 2j4d vs 1np7", "2j4d_1np7.png"),
    (["4dja.cif", "6k8i.cif"], "Degree Distribution: 4dja vs 6k8i", "4dja_6k8i.png"),
    (["1dnp.cif", "1iqr.cif"], "Degree Distribution: 1dnp vs 1iqr", "1dnp_1iqr.png")
]

# Plot and save figures
for protein_files, title, save_name in protein_sets:
    plot_degree_distributions(protein_files, title, save_name)


for cif_file in cif_files:
   _network(cif_file)
   _degree_distribution(cif_file)
