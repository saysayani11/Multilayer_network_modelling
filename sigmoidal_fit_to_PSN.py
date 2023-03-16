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
os.chdir("/mnt/808e843c-fb08-404a-9f45-efd817d81a77/Sayantoni/networks/cif_files_1758")

def read_textfile():
    with open('excluded','r') as text:
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

    prefix = "CA"
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
    return cutoff_list

def dict_pos(file):
    temp = fetch(file)
    res = temp[temp.columns[4]]
    xyz_matrix = temp[temp.columns[1:4]].astype(float)
   
    prefix = "CA"
    labels = [prefix + item for item in res]
    
    tuples =[]
    for item in combinations(labels,2):
        tuples.append(item)
       
    keys = labels
    dict_xyz = (xyz_matrix.astype(float).values.tolist())
    dict_network = dict(list((zip(keys, dict_xyz))))

    return (dict_network)

def percolation(file):
    
    dict_network = dict_pos(file)

    
    #-------------- PERCOLATION -------------
    
    indices = _network(file).index.values
    H = nx.Graph()
    largest_cc = []
    for i in range(1000):
        df_len = len(_network(file))      
        remove_frac = int((i/1000)*df_len)
        drop_indices = np.random.choice((indices), remove_frac, replace=False)
        df_subset  = _network(file).drop(drop_indices) 
        H= nx.from_pandas_edgelist(df_subset)    
        largest_cc.append(len(max(nx.connected_components(H), key=len))) 

        
    
    os.chdir("/mnt/808e843c-fb08-404a-9f45-efd817d81a77/Sayantoni/networks/cif_files_1758")
    #---Generate Percolation plot
    raw_y_data = list(largest_cc)
    ydata= (raw_y_data - (np.min(raw_y_data))) / (np.max(raw_y_data) - np.min(raw_y_data))
    xdata = np.arange(0, 1, 0.001).tolist()

    df = pd.DataFrame({'frac_of_removed_links': xdata,
                          'size_of_largest_component_normalized': ydata
                        }) 
    fig = plt.figure(figsize=(10,10))
    plt.plot(xdata,ydata, color='blue', linestyle='dashed', linewidth = 0.2,
          marker='o', markersize=3, label = 'Bond Percolation on PSN')
    
    plt.title(str(file) [:-4])
    plt.xlabel("fraction of removed links")
    plt.ylabel("size of the largest component (normalized)")
  
    
    #------------------------------------------ CURVE FITTING ----------------------------------
    from scipy.optimize import curve_fit
    def sigmoid(x, L, k, x0):
        y = L / (1 + np.exp(-k*(x-x0)))
        return y
    
    # Fit sigmoidal curve
    popt_sig, pcov_sig = curve_fit(sigmoid, xdata, ydata, method='lm', maxfev=100000)
    
    # Plot data and curves

    plt.plot(xdata, sigmoid(xdata, *popt_sig), 'r-', label='Sigmoidal Fit')
    plt.legend()
    
    #--------------------------------------------------------------------------------------------
    
    #---Save peroclation plot
    filename = str(file) [:-4]
    directory = "/mnt/808e843c-fb08-404a-9f45-efd817d81a77/Sayantoni/networks/results/" + filename
    subdirectory = directory + "/bond_percolation_100/percolation_plot"
    os.chdir(subdirectory)
    plt.savefig(filename + "_sigmoidfit" + ".tiff",dpi=300)
    os.chdir("/mnt/808e843c-fb08-404a-9f45-efd817d81a77/Sayantoni/networks/cif_files_1758")
    
    #---Save percolation plot data
    filename = str(file) [:-4]   
    directory = "/mnt/808e843c-fb08-404a-9f45-efd817d81a77/Sayantoni/networks/results/" + filename
    subdirectory = directory + "/bond_percolation_100/percolation_plot"
    os.chdir(subdirectory)
    df.to_excel(filename+ "_sigmoidfit" + ".xlsx")
    os.chdir("/mnt/808e843c-fb08-404a-9f45-efd817d81a77/Sayantoni/networks/cif_files_1758")
    

for i in cif_files:
    _network(i)
    percolation(i)
    
    