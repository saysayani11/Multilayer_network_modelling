import os
import numpy as np
import pandas as pd
from scipy import stats
from scipy.integrate import simps
from scipy.stats import gaussian_kde
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib import cm
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.MMCIFParser import MMCIFParser
from scipy.spatial.distance import pdist, squareform
import itertools
from itertools import combinations
from collections import Counter
import seaborn as sns
plt.rcParams.update({'figure.max_open_warning': 0, 'figure.dpi': 300})

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

def protein_data_within_range(file, residue_ranges):
    
   #--Count nodes and edges within given residue ranges for a protein structure
    node_data = fetch(file)
    edge_data = _network(file)
    
    node_counts = []
    edge_counts = []
    
    for start, stop in residue_ranges:
        nodes_in_range = node_data[(node_data['residue_no'].astype(int) >= start) & (node_data['residue_no'].astype(int) <= stop)]
        node_counts.append(len(nodes_in_range))
        
        edges_in_range = edge_data[edge_data['source'].str.extract('(\d+)')[0].astype(int).between(start, stop) & 
                                   edge_data['target'].str.extract('(\d+)')[0].astype(int).between(start, stop)]
        edge_counts.append(len(edges_in_range))
    
    result_df = pd.DataFrame({
        'Residue_Range_Start': [i[0] for i in residue_ranges],
        'Residue_Range_Stop': [i[1] for i in residue_ranges],
        'Node_Count': node_counts,
        'Edge_Count': edge_counts
    })
    
    return result_df

def consolidated_data_within_range(cif_files, ABD_ranges, FADBD_ranges):
    results = []
    
    for idx, file in enumerate(cif_files):
        ABD_result = protein_data_within_range(file, [ABD_ranges[idx]])
        ABD_result["File"] = os.path.basename(file)
        ABD_result["Domain"] = "ABD"
        
        FADBD_result = protein_data_within_range(file, [FADBD_ranges[idx]])
        FADBD_result["File"] = os.path.basename(file)
        FADBD_result["Domain"] = "FAD-BD"
        
        results.append(ABD_result)
        results.append(FADBD_result)
    
    consolidated_df = pd.concat(results, ignore_index=True)
    consolidated_df = consolidated_df[['File', 'Domain', 'Residue_Range_Start', 'Residue_Range_Stop', 'Node_Count', 'Edge_Count']]
    
    return consolidated_df

ABD_ranges = [ (5, 223), (5, 227), (38, 207), (19, 180), (25, 182), (6, 182), (5, 172), (6, 168), (3, 166), (3, 162), (21, 175), (6, 168), (179, 353), (84, 250), (3, 164), (2, 168), (6, 165), (12, 169), (5, 160) ]
FADBD_ranges = [ (226, 491), (229, 488), (234, 488), (228, 461), (212, 459), (213, 514), (290, 487), (297, 494), (286, 483), (288, 486), (306, 504), (288, 478), (393, 508), (368, 526), (241, 395), (272, 466), (274, 473), (290, 488), (287, 485) ]

cif_files = list(read_textfile())
consolidated_results = consolidated_data_within_range(cif_files, ABD_ranges, FADBD_ranges)
ABD_data = consolidated_results[consolidated_results['Domain'] == 'ABD']
FADBD_data = consolidated_results[consolidated_results['Domain'] == 'FAD-BD']

def calculate_closeness(file):
    edge_data = _network(file)
    node_data = fetch(file)

    # Create the graph
    G = nx.Graph()

    for _, row in node_data.iterrows():
        G.add_node(row['ca_atom_number'], label=row['residue_no'], residue_name=row['residue_name'])
        
    for _, row in edge_data.iterrows():
        G.add_edge(row['source'].split("ATM")[1], row['target'].split("ATM")[1], weight=row['distance'])


    # Calculate node closeness
    node_closeness = nx.closeness_centrality(G)

    # Generate node dataframe
    node_df_data = []
    for node, closeness in node_closeness.items():
        residue_id = G.nodes[node]['label']
        residue_name = G.nodes[node]['residue_name']
        node_df_data.append([residue_id, residue_name, closeness])

    node_df = pd.DataFrame(node_df_data, columns=['residue_id', 'residue_name', 'node_closeness'])
    node_df = node_df.sort_values(by='node_closeness', ascending=False)  # Sorting


    return node_df


def plot_closeness_distribution(cif_files):
    sns.set(style="white")
    sns.set_palette("flare")
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    for idx, file in enumerate(cif_files):
        protein_name = os.path.basename(file).split('.')[0]
        
        node_df = calculate_closeness(file)  # Corrected this line
        counts, bin_edges = np.histogram(node_df['node_closeness'], bins=20)
        ax.scatter(bin_edges[:-1], counts, label=protein_name, s=30)
    
    ax.set_xlabel("Closeness Centrality")
    ax.set_ylabel("Frequency")
    ax.set_title("Closeness Centrality Distribution")
    ax.legend()
    ax.grid(False)  # Turn off grids
    
    plt.tight_layout()
    plt.savefig("closeness_distribution.png", dpi=300)  # Save with 300 DPI
    plt.show()



# Call the function to plot the distribution
plot_closeness_distribution(cif_files)
