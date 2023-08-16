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

ABD_ranges = [ (5, 223), (5, 227), (38, 207), (19, 180), (31, 161), (6, 182), (5, 172), (6, 168), (3, 166), (3, 162), (21, 175), (6, 168), (179, 353), (84, 250), (3, 164), (2, 168), (6, 165), (12, 169), (5, 160) ]
FADBD_ranges = [ (226, 491), (229, 488), (234, 488), (228, 461), (252, 436), (213, 514), (290, 487), (297, 494), (286, 483), (288, 486), (306, 504), (288, 478), (393, 508), (368, 526), (241, 395), (272, 466), (274, 473), (290, 488), (287, 485) ]

cif_files = list(read_textfile())
consolidated_results = consolidated_data_within_range(cif_files, ABD_ranges, FADBD_ranges)
ABD_data = consolidated_results[consolidated_results['Domain'] == 'ABD']
FADBD_data = consolidated_results[consolidated_results['Domain'] == 'FAD-BD']




def plot_domain_network(file, ABD_range, FADBD_range):
    #-- Generate the complete network
    edge_data = _network(file)
    node_data = fetch(file)
    
   
    G = nx.Graph()
    
    #-- Add nodes and edges to the graph
    for _, row in node_data.iterrows():
        G.add_node(row['ca_atom_number'], label=row['residue_no'])
        
    for _, row in edge_data.iterrows():
        G.add_edge(row['source'].split("ATM")[1], row['target'].split("ATM")[1], weight=row['distance'])
    
    #-- Set the colors based on the domain
    node_colors = []
    edge_colors = []
    
    for node in G.nodes():
        residue_no = G.nodes[node]['label']
        
        if int(residue_no) >= ABD_range[0] and int(residue_no) <= ABD_range[1]:
            node_colors.append('red')
        elif int(residue_no) >= FADBD_range[0] and int(residue_no) <= FADBD_range[1]:
            node_colors.append('darkblue')
        else:
            node_colors.append('grey')  # default color
    
    for edge in G.edges():
        source_residue_no = G.nodes[edge[0]]['label']
        target_residue_no = G.nodes[edge[1]]['label']
        
        if int(source_residue_no) >= ABD_range[0] and int(source_residue_no) <= ABD_range[1] and int(target_residue_no) >= ABD_range[0] and int(target_residue_no) <= ABD_range[1]:
            edge_colors.append('pink')
        elif int(source_residue_no) >= FADBD_range[0] and int(source_residue_no) <= FADBD_range[1] and int(target_residue_no) >= FADBD_range[0] and int(target_residue_no) <= FADBD_range[1]:
            edge_colors.append('blue')
        else:
            edge_colors.append('lightgrey')  # default color
            
  


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


def _degree_with_domain(file, ABD_range, FADBD_range):
    edge_data = _network(file)
    node_data = fetch(file)

    # Create the graph
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
        
        # Determine if the residue belongs to ABD or FAD-BD
        if ABD_range[0] <= int(residue_id) <= ABD_range[1]:
            domain = 'ABD'
        elif FADBD_range[0] <= int(residue_id) <= FADBD_range[1]:
            domain = 'FAD-BD'
        else:
            domain = 'Other'
            
        df_data.append([residue_id, residue_name, degree, domain])

    df = pd.DataFrame(df_data, columns=['residue_id', 'residue_name', 'degree', 'domain'])
    
    return df

#-- List to store the dataframes
degree_dfs = []

#-- Loop over each CIF file
for file_path, ABD_range, FADBD_range in zip(cif_files, ABD_ranges, FADBD_ranges):
    # Compute degree dataframe with domain information
    degree_df = _degree_with_domain(file_path, ABD_range, FADBD_range)
    
    #-- Append the dataframe to the list
    degree_dfs.append(degree_df)

#-- 'degree_dfs' contains 18 dataframes, each for a different protein

def separate_dataframes_by_domain(degree_dfs):
    domain_data = {'ABD': [], 'FAD-BD': [], 'Other': []}
    
    for df in degree_dfs:
        for domain_type in ['ABD', 'FAD-BD', 'Other']:
            filtered_df = df[df['domain'] == domain_type]
            domain_data[domain_type].append(filtered_df)
    
    #-- Convert lists of dataframes to dictionaries
    for domain_type in domain_data:
        domain_data[domain_type] = {idx: df for idx, df in enumerate(domain_data[domain_type])}
    
    return domain_data

#-- Separate the dataframes by domain
separated_data = separate_dataframes_by_domain(degree_dfs)

# 'separated_data' is a dictionary with keys 'ABD', 'FAD-BD', and 'Other',
# where each key points to a dictionary of dataframes for that domain type

def apply_z_score_cutoff(dataframes_dict, z_cutoff=2):
    filtered_data = {}

    for domain_type, dfs_dict in dataframes_dict.items():
        filtered_data[domain_type] = {}

        for idx, df in dfs_dict.items():
            z_scores = stats.zscore(df['degree'])
            filtered_df = df[abs(z_scores) <= z_cutoff].copy()
            filtered_df['z_score'] = z_scores  # Add Z score column to the dataframe
            filtered_data[domain_type][idx] = filtered_df

    return filtered_data

# Apply Z score cutoff on the separated dataframes and include Z scores
filtered_data_with_z_scores = apply_z_score_cutoff(separated_data)

# 'filtered_data_with_z_scores' is a dictionary with the same structure as 'separated_data',
# but the dataframes within each domain type are filtered and include Z scores


def write_dataframes_to_excel(dataframes_dict, output_folder):
    os.makedirs(output_folder, exist_ok=True)
    
    for domain_type, dfs_dict in dataframes_dict.items():
        domain_folder = os.path.join(output_folder, domain_type)
        os.makedirs(domain_folder, exist_ok=True)
        
        for idx, df in dfs_dict.items():
            excel_path = os.path.join(domain_folder, f'dataframe_{domain_type}_{idx}.xlsx')
            with pd.ExcelWriter(excel_path) as writer:
                df.to_excel(writer, sheet_name=f'dataframe_{idx}', index=False)

# Provide the path to the output folder where Excel files will be saved
output_folder = r"C:\Users\saySa\OneDrive\Desktop\task_1\PSNs\networks"

# Write dataframes with Z scores to separate Excel files
write_dataframes_to_excel(filtered_data_with_z_scores,output_folder)

