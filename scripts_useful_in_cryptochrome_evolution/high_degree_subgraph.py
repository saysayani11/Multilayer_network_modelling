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

def _degree(file):
    edge_data = _network(file)
    node_data = fetch(file)

    #-- Create the graph
    G = nx.Graph()

    for _, row in node_data.iterrows():
        G.add_node(row['ca_atom_number'], label=row['residue_no'], residue_name=row['residue_name'])
        
    for _, row in edge_data.iterrows():
        G.add_edge(row['source'].split("ATM")[1], row['target'].split("ATM")[1], weight=row['distance'])

    #-- Calculate degree for each node
    node_degrees = G.degree()

    #-- Generate dataframe
    df_data = []
    for node, degree in node_degrees:
        residue_id = G.nodes[node]['label']
        residue_name = G.nodes[node]['residue_name']
        df_data.append([residue_id, residue_name, degree])

    df = pd.DataFrame(df_data, columns=['residue_id', 'residue_name', 'degree'])
       
    return df


def create_high_degree_subgraphs(file):
    # Calculate node degrees
    df_degree = _degree(file)
    
    # Calculate the 90th percentile of the degree distribution
    degree_threshold = df_degree['degree'].quantile(0.9)
    
    # Select nodes with degree higher than the 90th percentile
    high_degree_nodes = df_degree[df_degree['degree'] > degree_threshold]
    # print("High-degree nodes:")
    # print(high_degree_nodes)
    
    # # Extract the protein name or ID from the file variable
    # protein_name = file.split('.')[0]  # assuming the file name is 'proteinID.cif'
    
    # # Create a unique file name for each protein
    # output_filename = f"{protein_name}_high_degree_nodes.xlsx"
    
    # # Save the high_degree_nodes DataFrame to an Excel file
    # high_degree_nodes.to_excel(output_filename, index=False)
    
    return high_degree_nodes



def _subgraph(file):
    # Calculate the high-degree nodes
    high_degree_nodes = create_high_degree_subgraphs(file)
    
    # Get the original graph
    edge_data = _network(file)
    
    # Fetch the complete atom data using the fetch function
    atom_data = fetch(file)
    
    # Filter the atom data to retain only those atoms corresponding to high-degree residues
    high_degree_residues = high_degree_nodes['residue_id'].astype(str).tolist()
    high_degree_atoms = atom_data[atom_data['residue_no'].isin(high_degree_residues)]
    
    # Extract the atom_numbers corresponding to the high-degree atoms
    high_degree_atom_numbers = high_degree_atoms['ca_atom_number'].astype(str).tolist()
    
    # Prefix the atom numbers with "ATM" to match the source and target columns format in edge_data
    high_degree_atom_ids = ["ATM" + num for num in high_degree_atom_numbers]
    
    # Filter the original edge data to keep only those rows where both source and target are in high_degree_atom_ids
    subgraph_edge_data = edge_data[edge_data['source'].isin(high_degree_atom_ids) & edge_data['target'].isin(high_degree_atom_ids)]
    
    # Mapping the source and target back to residue_id
    atom_to_residue = high_degree_atoms.set_index('ca_atom_number')['residue_no'].to_dict()
    subgraph_edge_data['source_residue_id'] = subgraph_edge_data['source'].str.extract('(\d+)')[0].map(atom_to_residue)
    subgraph_edge_data['target_residue_id'] = subgraph_edge_data['target'].str.extract('(\d+)')[0].map(atom_to_residue)
    
    # Reordering the columns
    subgraph_edge_data = subgraph_edge_data[['source_residue_id', 'target_residue_id', 'distance']]
    
    # Extract the protein name or ID from the file variable
    protein_name = file.split('.')[0]  # assuming the file name is 'proteinID.cif'
    
    # Create a unique file name for the subgraph and node information
    output_filename_subgraph = f"{protein_name}_subgraph.xlsx"
    output_filename_nodes = f"{protein_name}_subgraph_nodes.xlsx"
    
    # Save the modified subgraph edge data DataFrame to an Excel file
    subgraph_edge_data.to_excel(output_filename_subgraph, index=False)
    
    # For each node in the subgraph, get the distance information and save it
    subgraph_nodes = high_degree_atoms[['residue_no', 'residue_name']]
    subgraph_nodes['distance'] = subgraph_edge_data.groupby('source_residue_id')['distance'].sum().reset_index()['distance']
    
    # Save the node information for the subgraph to an Excel file
    subgraph_nodes.to_excel(output_filename_nodes, index=False)
    
    return subgraph_edge_data, subgraph_nodes


# Loop through each CIF file and create high degree nodes subgraphs
for file in cif_files:
    _network(file)
    _degree(file)
    create_high_degree_subgraphs(file)
    _subgraph(file)
 
