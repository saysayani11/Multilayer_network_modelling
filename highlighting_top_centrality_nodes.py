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

# Set Path
os.chdir("/home/sayantoni/Desktop/")
file = "3zxs.cif"

# Read values from Excel file
values_file = "values.xlsx"  
values_data = pd.read_excel(values_file)

atom_type = MMCIF2Dict(file)['_atom_site.group_PDB']
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


# --- Pick CA
pick_CA = data_cif[data_cif[3] == 'CA']

# --- choose chain A
pick_CA_chainA = pick_CA[pick_CA[6] == 'A']

#--- Pick XYZ coordinates
xyz_coord = pick_CA_chainA[pick_CA_chainA.columns[9:12]]
   
#--- Pick the residue numbers
res =  pick_CA_chainA[7]

#--- Pick atomsite IDs  
atomid = pick_CA_chainA[1]

# --- AssymID
assymid = pick_CA_chainA[6]

# --- Residue name
resname = pick_CA_chainA[5]

temp = pd.concat([res, xyz_coord, atomid, assymid, resname], axis=1)
temp2 = temp


res = temp[temp.columns[4]]
xyz_matrix = temp[temp.columns[1:4]]
points = np.array(xyz_matrix).astype(float)
dist_condensed = pdist(points)

prefix = "ATM"
labels = [prefix + item for item in res]

tuples = []
for item in combinations(labels, 2):
    tuples.append(item)

source = pd.DataFrame(tuples)[0]
target = pd.DataFrame(tuples)[1]

edgelist = pd.DataFrame({'source': source,
                          'target': target,
                          'distance': dist_condensed
                          })

filename = str(file)
cutoff_list = edgelist[edgelist['distance']<8]

# --- SET COORDINATES FOR THE NETWORK

chains =  list(set(temp[temp.columns[5]] ))
count_chain_res = Counter(temp[temp.columns[5]])

keys = labels
dict_xyz = (xyz_matrix.astype(float).values.tolist())
dict_network = dict(list((zip(keys, dict_xyz))))


# --- PLOT THE NETWORK
H = nx.Graph()
H = nx.from_pandas_edgelist(cutoff_list)
pos = dict_network
temp = list(H.edges())

edge_xyz = [(pos[u], pos[v]) for u, v in H.edges()]
t = np.asarray(edge_xyz)

fig = plt.figure(figsize=(20,20))
ax = fig.add_subplot(111, projection="3d")
plt.rcParams['figure.dpi'] = 300
   
# xs = np.array(xyz_matrix[9],dtype=float)
# ys = np.array(xyz_matrix[10],dtype=float)
# zs = np.array(xyz_matrix[11],dtype=float)
# ax.scatter(xs,ys,zs,color='pink')



#-------------highlighting---------------
subset = temp2[temp2.columns[1:5]].astype(float)
kk = list(values_data["ca_atom"])               
selected_df = subset[subset.iloc[:, 3].isin(kk)]
            
xs_subset = np.array(selected_df[9],dtype=float)
ys_subset = np.array(selected_df[10],dtype=float)
zs_subset = np.array(selected_df[11],dtype=float)
ax.scatter(xs_subset,ys_subset,zs_subset,color='red')

cutoff_list['source'] = cutoff_list['source'].str[3:].astype(int)
cutoff_list['target'] = cutoff_list['target'].str[3:].astype(int)
mask = cutoff_list['source'].isin(kk) | cutoff_list['target'].isin(kk)
df_filtered = cutoff_list[mask]

H = nx.Graph()
H = nx.from_pandas_edgelist(cutoff_list)

edge_xyz = [(pos[u], pos[v]) for u, v in H.edges()]
t = np.asarray(edge_xyz)

ax.grid(True)
ax.set_axis_off()
for vizedge in t:
    ax.plot(*vizedge.T, linewidth=0.8,color="green")
plt.title(str(file))  
filename = str(file)




