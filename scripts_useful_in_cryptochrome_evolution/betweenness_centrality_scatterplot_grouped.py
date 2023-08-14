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
            
    #-- Plot the network
    fig, ax = plt.subplots(figsize=(10, 10))
    pos = nx.spring_layout(G)
    nx.draw(G, pos, ax=ax, with_labels=False, node_size=100, node_color=node_colors, edge_color=edge_colors, width=0.5)
    ax.set_title(os.path.basename(file))
    plt.show()


    # output_filename = f"{os.path.basename(file).split('.')[0]}_network.tiff"
    # fig.savefig(output_filename, dpi=300, format='tiff')
    # print(f"Saved plot as: {output_filename}")
    # plt.close(fig)
    
############################ PLOT ##################################

#-- Plotting for each protein structure manually
plot_domain_network(cif_files[0], ABD_ranges[0], FADBD_ranges[0])
plot_domain_network(cif_files[1], ABD_ranges[1], FADBD_ranges[1])
plot_domain_network(cif_files[2], ABD_ranges[2], FADBD_ranges[2])
plot_domain_network(cif_files[3], ABD_ranges[3], FADBD_ranges[3])
plot_domain_network(cif_files[4], ABD_ranges[4], FADBD_ranges[4])
plot_domain_network(cif_files[5], ABD_ranges[5], FADBD_ranges[5])
plot_domain_network(cif_files[6], ABD_ranges[6], FADBD_ranges[6])
plot_domain_network(cif_files[7], ABD_ranges[7], FADBD_ranges[7])
plot_domain_network(cif_files[8], ABD_ranges[8], FADBD_ranges[8])
plot_domain_network(cif_files[9], ABD_ranges[9], FADBD_ranges[9])
plot_domain_network(cif_files[10], ABD_ranges[10], FADBD_ranges[10])
plot_domain_network(cif_files[11], ABD_ranges[11], FADBD_ranges[11])
plot_domain_network(cif_files[12], ABD_ranges[12], FADBD_ranges[12])
plot_domain_network(cif_files[13], ABD_ranges[13], FADBD_ranges[13])
plot_domain_network(cif_files[14], ABD_ranges[14], FADBD_ranges[14])
plot_domain_network(cif_files[15], ABD_ranges[15], FADBD_ranges[15])
plot_domain_network(cif_files[16], ABD_ranges[16], FADBD_ranges[16])
plot_domain_network(cif_files[17], ABD_ranges[17], FADBD_ranges[17])

######################################################################

############# CALCULATE BETWENNESS CENTRALITY  #######################
def calculate_betweenness(file):
    edge_data = _network(file)
    node_data = fetch(file)

    # Create the graph
    G = nx.Graph()

    for _, row in node_data.iterrows():
        G.add_node(row['ca_atom_number'], label=row['residue_no'], residue_name=row['residue_name'])
        
    for _, row in edge_data.iterrows():
        G.add_edge(row['source'].split("ATM")[1], row['target'].split("ATM")[1], weight=row['distance'])

    # Calculate edge betweenness
    edge_betweenness = nx.edge_betweenness_centrality(G, normalized=True)

    # Calculate node betweenness
    node_betweenness = nx.betweenness_centrality(G, normalized=True)

    # Generate node dataframe
    node_df_data = []
    for node, betweenness in node_betweenness.items():
        residue_id = G.nodes[node]['label']
        residue_name = G.nodes[node]['residue_name']
        node_df_data.append([residue_id, residue_name, betweenness])

    node_df = pd.DataFrame(node_df_data, columns=['residue_id', 'residue_name', 'node_betweenness'])
    node_df = node_df.sort_values(by='node_betweenness', ascending=False)  # Sorting

    # Generate edge dataframe
    edge_df_data = []
    for edge, betweenness in edge_betweenness.items():
        source_residue_id = G.nodes[edge[0]]['label']
        target_residue_id = G.nodes[edge[1]]['label']
        edge_df_data.append([f"{source_residue_id}-{target_residue_id}", betweenness])

    edge_df = pd.DataFrame(edge_df_data, columns=['residue_pair', 'edge_betweenness'])
    edge_df = edge_df.sort_values(by='edge_betweenness', ascending=False)  # Sorting

    return node_df, edge_df

# Modified saving loop to save sorted dataframes
# for file in cif_files:
#     protein_name = os.path.basename(file).split('.')[0]
    
#     node_df, edge_df = calculate_betweenness(file)
    
#     node_save_path = f"NodeBetweenness_{protein_name}_Sorted.csv"
#     edge_save_path = f"EdgeBetweenness_{protein_name}_Sorted.csv"
    
#     node_df.to_csv(node_save_path, index=False)
#     edge_df.to_csv(edge_save_path, index=False)

#     print(f"Saved {node_save_path} and {edge_save_path}!")


######################################################################


def _degree(file):
    edge_data = _network(file)
    node_data = fetch(file)

    # Create the graph
    G = nx.Graph()

    for _, row in node_data.iterrows():
        G.add_node(row['ca_atom_number'], label=row['residue_no'], residue_name=row['residue_name'])
        
    for _, row in edge_data.iterrows():
        G.add_edge(row['source'].split("ATM")[1], row['target'].split("ATM")[1], weight=row['distance'])

    # Calculate degree for each node
    node_degrees = G.degree()

    # Generate dataframe
    df_data = []
    for node, degree in node_degrees:
        residue_id = G.nodes[node]['label']
        residue_name = G.nodes[node]['residue_name']
        df_data.append([residue_id, residue_name, degree])

    df = pd.DataFrame(df_data, columns=['residue_id', 'residue_name', 'degree'])
    
    # Calculate Z-scores
    mean_degree = df['degree'].mean()
    std_degree = df['degree'].std()
    df['z_score'] = (df['degree'] - mean_degree) / std_degree
    
    # Applying the cutoffs
    df_2 = df[df['z_score'].abs() >= 2]
    df_3 = df[df['z_score'].abs() >= 3]
    
    return df, df_2, df_3


# # Saving the calculated and sorted degrees for each protein structure
# for file in cif_files:
#     protein_name = os.path.basename(file).split('.')[0]
    
#     df, df_2, df_3 = _degree(file)
    
#     save_path = f"NodeDegree_{protein_name}.csv"
#     df.to_csv(save_path, index=False)

#     save_path_2 = f"NodeDegree_ZScore2_{protein_name}.csv"
#     df_2.to_csv(save_path_2, index=False)

#     save_path_3 = f"NodeDegree_ZScore3_{protein_name}.csv"
#     df_3.to_csv(save_path_3, index=False)

#     print(f"Saved {save_path}, {save_path_2}, and {save_path_3}!")
    
def plot_betweenness_distribution(cif_files):
    sns.set(style="white")
    sns.set_palette("flare")
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    for idx, file in enumerate(cif_files):
        protein_name = os.path.basename(file).split('.')[0]
        
        node_df, _ = calculate_betweenness(file)
        counts, bin_edges = np.histogram(node_df['node_betweenness'], bins=20)
        ax.scatter(bin_edges[:-1], counts, label=protein_name, s=30)
    
    ax.set_xlabel("Betweenness Centrality")
    ax.set_ylabel("Frequency")
    ax.set_title("Betweenness Centrality Distribution")
    ax.legend()
    ax.grid(False)  # Turn off grids
    
    plt.tight_layout()
    plt.savefig("betweenness_distribution.png", dpi=300)  # Save with 300 DPI
    plt.show()

# Call the function to plot the distribution
plot_betweenness_distribution(cif_files)

def calculate_min_max_betweenness(cif_files):
    min_max_data = []
    
    for file in cif_files:
        protein_name = os.path.basename(file).split('.')[0]
        node_df, _ = calculate_betweenness(file)
        
        min_bc = node_df['node_betweenness'].min()
        max_bc = node_df['node_betweenness'].max()
        
        min_max_data.append([protein_name, min_bc, max_bc])
    
    min_max_df = pd.DataFrame(min_max_data, columns=['Protein', 'Min_Betweenness', 'Max_Betweenness'])
    return min_max_df

# Call the function to calculate and store the min-max data
min_max_betweenness_df = calculate_min_max_betweenness(cif_files)
print(min_max_betweenness_df)


def plot_min_max_betweenness(min_max_betweenness_df):
    sns.set_palette("tab10")
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    for _, row in min_max_betweenness_df.iterrows():
        ax.plot([row['Protein'], row['Protein']], [row['Min_Betweenness'], row['Max_Betweenness']], marker='o')
        ax.scatter(row['Protein'], row['Min_Betweenness'], color='blue', marker='o', label='Min', s=50)
        ax.scatter(row['Protein'], row['Max_Betweenness'], color='red', marker='o', label='Max', s=50)
    
    ax.set_xlabel("Protein")
    ax.set_ylabel("Betweenness Centrality")
    ax.set_title("Min and Max Betweenness Centrality for Proteins")
   
    ax.grid(False)  # Turn off grids
    
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig("min_max_betweenness.png", dpi=300)  # Save with 300 DPI
    plt.show()

# Call the function to plot the min-max data
plot_min_max_betweenness(min_max_betweenness_df)

def plot_betweenness_distribution_groups_subplots(cif_files, groups, group_names):
    sns.set(style="white")
    sns.set_palette("flare")
    
    num_groups = len(groups)
    fig, axs = plt.subplots(nrows=num_groups, figsize=(10, 6 * num_groups), sharex=True)
    
    for idx, (group, group_name) in enumerate(zip(groups, group_names)):
        axs[idx].set_ylabel("Frequency")
        axs[idx].set_title(f"Betweenness Centrality Distribution - {group_name}")
        axs[idx].grid(False)  # Turn off grids
        
        for file in group:
            protein_name = os.path.basename(file).split('.')[0]
            node_df, _ = calculate_betweenness(file)
            counts, bin_edges = np.histogram(node_df['node_betweenness'], bins=20)
            axs[idx].scatter(bin_edges[:-1], counts, label=protein_name, s=30)
        
        axs[idx].legend()
        
    axs[-1].set_xlabel("Betweenness Centrality")
    plt.tight_layout()
    plt.savefig("betweenness_distribution_subplots.png", dpi=300)  # Save with 300 DPI
    plt.show()

# Define groups and group_names
group1 = ['3ZXS.cif', '4DJA.cif']
group2 = ['3UMV.cif', '2XRY.cif', 'bcpd2_cut.cif']
group3 = ['4GU5.cif', '5ZM0.cif', '7AYV.cif', '6PTZ.cif', '5T5X.cif', '7V8Y.cif']
group4 = ['1NP7.cif', '2J4D.cif']
group5 = ['1IQR.cif', '1DNP.cif', '4U63.cif', '1U3C.cif', '6K8I.cif']

groups = [group1, group2, group3, group4, group5]
group_names = ["Group 1 (6-4 like /FeS)", "Group 2 (CPD 2)", "Group 3 (3.6-4/Animal)", "Group 4 (DASH)", "Group 5 (5.CPD1/3)"]

# Call the function to plot the distribution by groups in subplots
plot_betweenness_distribution_groups_subplots(cif_files, groups, group_names)

