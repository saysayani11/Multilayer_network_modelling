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
plt.rcParams.update({'figure.max_open_warning': 0})

#-- Set Path
os.chdir(r"C:\Users\saySa\OneDrive\Desktop\test")

def read_textfile():
    with open('reps.txt','r') as text:
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
   
    
    #---  remove extra chains
    remove_chains = data_cif [data_cif[6] == 'A']
    
    #--- Pick CA
    pick_CA =  remove_chains[remove_chains[3] == 'CA']
   
    #--- Pick XYZ coordinates
    xyz_coord = pick_CA[pick_CA.columns[9:12]]
   
    #--- Pick the residue numbers
    res =  pick_CA[7]
   
    #--- Pick atomsite IDs  
    atomid = pick_CA[1]
   
    #--- AssymID
    assymid = pick_CA[6]
   
    #--- Residue name
    resname = pick_CA[5]
    temp = pd.concat([res, xyz_coord, atomid, assymid, resname], axis=1)
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
    cutoff_list = edgelist[edgelist['distance']<8]
     
    return cutoff_list

def _BCA(file):
    
    #--- Calculate normalized betweeness
    temp = _network(file)
    H = nx.Graph()
    H = nx.from_pandas_edgelist(temp)
    bca = nx.betweenness_centrality(H,normalized = True)  
    df = pd.DataFrame(list(bca.items()), columns=['Key', 'Value'])
    df['Key'] = df['Key'].str[3:].astype(int)
    df_sorted = df.sort_values(by='Key', ascending=True)
    
    
    #--- fetch and annotate the corresponding residue IDs
    parser = MMCIFParser()
    structure = parser.get_structure(str(file)[:4], str(file))
    model = structure[0]
    chain = model['A']
    residue_ids = []
    residue_ids = []
    for atom in chain.get_atoms():
        if atom.get_name() == 'CA':
            residue_id = atom.get_parent().get_id()[1]
            residue_ids.append(residue_id)
                
    #--- define a dataframe with (a) Residue ID (b) Atom ID and (c) Betweeness (normalized)              
    final_df = pd.concat([df_sorted, pd.DataFrame(residue_ids)], axis=1) 
    final_df = final_df.rename(columns={'Atom_ID': 'X', 'Value': 'Normalized_Betweeness' , 0 : 'Residue_ID'})
    # pd.DataFrame(residue_ids).to_excel("resID.xlsx")
    df_sorted.to_excel("2j4d.xlsx")
    return final_df
 
def histogram_BCA(file):
    y_data =  _BCA(file)['Normalized_Betweeness']
    plt.figure(dpi=300)
    n, bins, patches = plt.hist(y_data, bins=50, color='mediumvioletred')  
    for patch in patches:
        patch.set_edgecolor('white')
    plt.title('Betweenness Centrality Histogram')
    plt.xlabel('Betweenness Centrality Values')
    plt.ylabel('Frequency')
    plt.show() 
    return y_data

def KDE_BCA(file):
    y_data = _BCA(file)['Normalized_Betweeness']   
    plt.figure(dpi=300)
    sns.kdeplot(data=y_data, fill=True, common_norm=False, palette="crest", alpha=.5, linewidth=0)   
    plt.title('Betweenness Centrality KDE Plot')
    plt.xlabel('Betweenness Centrality Values')
    plt.ylabel('Density')  
    plt.show()    
    return y_data

def KDE_BCA():
    y_data_1 = _BCA("2j4d.cif")['Normalized_Betweeness']
    y_data_2 = _BCA("3umv.cif")['Normalized_Betweeness']
    y_data_3 = _BCA("3zxs.cif")['Normalized_Betweeness']
    y_data_4 = _BCA("4u63.cif")['Normalized_Betweeness']
    y_data_5 = _BCA("4gu5.cif")['Normalized_Betweeness']   
    fig, axes = plt.subplots(5, 1, figsize=(8, 20), dpi=300, sharex=True)
    sns.kdeplot(data=y_data_1, ax=axes[0], fill=True, common_norm=False, palette="crest", alpha=.5, linewidth=0, color='blue')
    sns.kdeplot(data=y_data_2, ax=axes[1], fill=True, common_norm=False, palette="crest", alpha=.5, linewidth=0, color='green')
    sns.kdeplot(data=y_data_3, ax=axes[2], fill=True, common_norm=False, palette="crest", alpha=.5, linewidth=0, color='red')
    sns.kdeplot(data=y_data_4, ax=axes[3], fill=True, common_norm=False, palette="crest", alpha=.5, linewidth=0, color='purple')
    sns.kdeplot(data=y_data_5, ax=axes[4], fill=True, common_norm=False, palette="crest", alpha=.5, linewidth=0, color='orange')

    legend_labels = ['2j4d', '3umv', '3zxs', '4u63', '4gu5']
    fig.legend(title='PSNs', labels=legend_labels, loc='upper center')
    
    fig.text(0.04, 0.5, 'Density', va='center', rotation='vertical')
    plt.savefig('kde_plots.tiff', dpi=300, format='tiff')
    plt.show() 
    return y_data_1, y_data_2, y_data_3, y_data_4, y_data_5

def histogram_BCA_line_plots():
    y_data_1 = _BCA("2j4d.cif")['Normalized_Betweeness']
    y_data_2 = _BCA("3umv.cif")['Normalized_Betweeness']
    y_data_3 = _BCA("3zxs.cif")['Normalized_Betweeness']
    y_data_4 = _BCA("4u63.cif")['Normalized_Betweeness']
    y_data_5 = _BCA("4gu5.cif")['Normalized_Betweeness'] 
    plt.figure(figsize=(10, 6), dpi=300)
    
    sns.kdeplot(data=y_data_1, fill=None, common_norm=False, palette="crest", alpha=.5, linewidth=1, color='blue')
    sns.kdeplot(data=y_data_2, fill=None, common_norm=False, palette="crest", alpha=.5, linewidth=1, color='green')
    sns.kdeplot(data=y_data_3, fill=None, common_norm=False, palette="crest", alpha=.5, linewidth=1, color='red')
    sns.kdeplot(data=y_data_4, fill=None, common_norm=False, palette="crest", alpha=.5, linewidth=1, color='purple')
    sns.kdeplot(data=y_data_5, fill=None, common_norm=False, palette="crest", alpha=.5, linewidth=1, color='orange')    
    plt.title('Betweenness Centrality KDE Plots', fontsize=16)
    plt.xlabel('Betweenness Centrality Values')
    plt.ylabel('Density')   
    legend_labels = ['2j4d', '3umv', '3zxs', '4u63', '4gu5']
    plt.legend(title='PSNs', labels=legend_labels)    
    plt.show()
    
    return y_data_1, y_data_2, y_data_3, y_data_4, y_data_5

def KDE_comparison():
    
#     #--- Compute peak values
#     y_data_1 = _BCA("4gu5.cif")['Normalized_Betweeness']   # Replace the PDB ID
#     kde_1 = sns.kdeplot(data=y_data_1)    
#     peak_index_1 = np.argmax(kde_1.get_lines()[0].get_data()[1])
#     peak_value_1 = kde_1.get_lines()[0].get_data()[0][peak_index_1]
#     print("COMPUTING PEAK VALUES:")
#     print(peak_value_1)
        
#     #--- Compute overlapping areas       
#     data = np.array([_BCA("2j4d.cif")['Normalized_Betweeness'],
#                       _BCA("3umv.cif")['Normalized_Betweeness'],
#                       _BCA("3zxs.cif")['Normalized_Betweeness'],
#                       _BCA("4u63.cif")['Normalized_Betweeness'],
#                       _BCA("4gu5.cif")['Normalized_Betweeness']],dtype=object)
#     data = [np.array(d)[np.isfinite(d)] for d in data]
#     overlaps = []
#     for i, j in itertools.combinations(range(len(data)), 2):
#         kde1 = gaussian_kde(data[i])
#         kde2 = gaussian_kde(data[j])
#         common_x_range = np.linspace(min(min(data[i]), min(data[j])), max(max(data[i]), max(data[j])), 1000)
#         y1 = kde1(common_x_range)
#         y2 = kde2(common_x_range)
#         overlapping_area = simps(np.minimum(y1, y2), common_x_range)
#         overlaps.append(overlapping_area)
#         print(f"Overlap between dataset {i+1} and {j+1}: {overlapping_area:.3f}")
#     return data
        
        






















    

 
 

