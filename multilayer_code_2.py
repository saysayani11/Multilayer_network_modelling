import os
import networkx as nx
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from scipy.spatial.distance import pdist, squareform
import pandas as pd
import numpy as np
from itertools import combinations
import matplotlib.pyplot as plt
from collections import Counter


#-- Set Path to PDB Files
os.chdir("/home/sayantoni/Desktop/test_pp_complex")

#-- READ MMCIF FILE

parser = MMCIFParser()
structure = parser.get_structure("1us7", "1us7.cif")
file_name = "1us7.cif"

atom_type = MMCIF2Dict(file_name) ['_atom_site.group_PDB']
atom_site = MMCIF2Dict(file_name) ['_atom_site.id']
atom_fullsymbol = MMCIF2Dict(file_name) ['_atom_site.label_atom_id']
residue_name = MMCIF2Dict(file_name) ['_atom_site.label_comp_id']
atom_assym = MMCIF2Dict(file_name) ['_atom_site.label_asym_id']
residue_id = MMCIF2Dict(file_name) ['_atom_site.label_seq_id']
entity_id = MMCIF2Dict(file_name) ['_atom_site.label_entity_id']
x_coord = MMCIF2Dict(file_name) ['_atom_site.Cartn_x']
y_coord = MMCIF2Dict(file_name) ['_atom_site.Cartn_y']
z_coord = MMCIF2Dict(file_name) ['_atom_site.Cartn_z']


#-- CREATE DATAFRAME
data_cif = pd.DataFrame(list(zip(atom_type, 
                    atom_site,
                    atom_fullsymbol,
                    residue_name, 
                    atom_assym,
                    residue_id,
                    entity_id,
                    x_coord, 
                    y_coord, 
                    z_coord)))

data_cif.columns = ['atom_type', 
                    'atom_site',
                    'atom_fullsymbol',
                    'residue_name', 
                    'atom_assym',
                    'residue_id',
                    'entity_id',
                    'x_coord', 
                    'y_coord', 
                    'z_coord']

    
#--- remove HETATM, TER
data_cif = (data_cif[~data_cif["atom_type"].isin(['HETATM','TER'])])

#--- PICK CA
pick_CA =  data_cif[data_cif["atom_fullsymbol"] == 'CA']

#--- FIND OUT NUMBER OF DISCINCT PROTEIN ENTITIES
distinct_entities = list(set(pick_CA["entity_id"]))

#--- CONVERT DATAFRAME INTO A EDGELIST
res = pick_CA[pick_CA.columns[5]]
xyz_matrix = pick_CA[pick_CA.columns[7:10]]
points = np.array(xyz_matrix).astype(float)
dist_condensed = pdist(points)
distances = squareform(dist_condensed)

#--- LABELS(by residue ID)
labels = []
pick_res = pick_CA[pick_CA.columns[1]]
for i in pick_res:
    labels.append(i + "_CA")
    
    
#--- DIVIDE pick_CA based on entity_id
chain_1_df = pick_CA[pick_CA['entity_id'] == "1"]
chain_2_df = pick_CA[pick_CA['entity_id'] == "2"]
chain = pick_CA["entity_id"]

t = pd.DataFrame(chain)

t['group'] = t['entity_id'].ne(t['entity_id'].shift()).cumsum()
df = t.groupby('group')
dfs = []
for name, data in df:
    dfs.append(data)

last_element = dfs[0][-1:]
nn = last_element.index.tolist()
    


#--- Divide Labels into two based on res_id
labels_1 = labels[0:len(chain_1_df)]
labels_2 = labels[len(chain_1_df):]


#--- CREATE TUPLES (SOURCE, TARGET) BASED ON THE LABELS
tuples =[]
for i in combinations(labels,2):
        tuples.append(i)
        
tuples_1 = []
for i in combinations(labels_1,2):
    tuples_1.append(i)
    

        
#--- GENERATE EDGELIST
edgelist = pd.DataFrame({'source': pd.DataFrame(tuples)[0],
                              'target': pd.DataFrame(tuples)[1],
                              'weights': dist_condensed
                            })

edgelist_1 = edgelist[0:len(tuples_1)]
edgelist_2 = edgelist[len(tuples_1):]

#--- SET CUTOFF AT 8 ANGSTROM
cutoff_list = edgelist[edgelist['weights']<8]

cutofflist_1 = edgelist_1[edgelist_1['weights']<8]
cutofflist_2 = edgelist_2[edgelist_2['weights']<8]


breakpoint_numbers_source= []
for i in list(cutoff_list["source"]):
    breakpoint_numbers_source.append(i[:len(i) - 3])
    
breakpoint_numbers_target= []
for i in list(cutoff_list["target"]):
    breakpoint_numbers_target.append(i[:len(i) - 3])
    
new_edgelist = pd.DataFrame({'source'  : breakpoint_numbers_source,
                             'target'  : breakpoint_numbers_target}).astype(int)


df1 = new_edgelist [new_edgelist ['source'] <=nn[0]]
df1['source'] = df1['source'].astype(str)+'_CA'
df1['target'] = df1['target'].astype(str)+'_CA'
df2 = new_edgelist [new_edgelist ['source'] > nn[0]]
df2['source'] = df2['source'].astype(str)+'_CA'
df2['target'] = df2['target'].astype(str)+'_CA'




# data1 = cutoff_list["source"]
# data2 = cutoff_list["target"]

# values1 = labels_1
# values2 = labels_2



    
    #  for j in i:
    # (any(i in x  for x in labels_1) == True & any(i in x  for x in labels_2) == True):
    #     new_cutoff.append(i)




#--- GENERATE GRAPH OBJECT
H = nx.Graph()
H = nx.from_pandas_edgelist(cutoff_list)

H1 = nx.from_pandas_edgelist(df1)
H2 = nx.from_pandas_edgelist(df2)



#--- VISUALIZE VIA MATPLOTLIB
    #--- NODES SCATTER-PLOT
fig = plt.figure(figsize=(20,20))
ax = fig.add_subplot(111, projection="3d")
xs = np.array(xyz_matrix["x_coord"],dtype=float)
ys = np.array(xyz_matrix["y_coord"],dtype=float)
zs = np.array(xyz_matrix["z_coord"],dtype=float)
ax.scatter(xs,ys,zs)




dict_xyz = (xyz_matrix.astype(float).values.tolist())

dict_network = dict(list((zip(labels, dict_xyz))))



edge_xyz_1 = [(dict_network[u], dict_network[v]) for u, v in H1.edges()]
edge_xyz_2 = [(dict_network[u], dict_network[v]) for u, v in H2.edges()]



t1 = np.asarray(edge_xyz_1)
t2 = np.asarray(edge_xyz_2)



for vizedge in t1:
    ax.plot(*vizedge.T, linewidth=0.8,color="red")
    
for vizedge in t2:
    ax.plot(*vizedge.T, linewidth=0.8,color="green")    




 
