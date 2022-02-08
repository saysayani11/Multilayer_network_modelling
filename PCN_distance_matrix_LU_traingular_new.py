from Bio.PDB import PDBList
from Bio.PDB import MMCIFParser
from itertools import combinations
from matplotlib.pyplot import figure 
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import networkx as nx

#-- Read mmCIF Files
get_pdb = PDBList()
structure_input = input("enter PDB ID: ")
get_pdb.retrieve_pdb_file(structure_input, pdir = '.', file_format = 'mmCif')
parser = MMCIFParser(QUIET = True)
structure = parser.get_structure(structure_input, structure_input+".cif")
first_model = structure[0]
#-- Calculate 3D distances between atom-pairs
def pairwise_atom_distance(dist_list = [], *args):
    x = list(dist_list)
    x1 = x[0]
    y1 = x[1]
    z1 = x[2]
    x2 = x[3]
    y2 = x[4]
    z2 = x[5]
    import math
    distance = math.sqrt(((x2-x1)**2) + ((y2-y1)**2) + ((z2-z1)**2))
    return(distance)
    
#-- Count CA atoms
count = 0    

for chain in first_model:
    for residue in chain:
        if residue.has_id("CA"):
            count = count + 1
        else:
             break
           
            
#-- Return 3D Coordinates (ca_coord_array)
ca_coord=[]
for chain in first_model:
    for residue in chain:
        if residue.has_id("CA"):
           ca = residue["CA"]
           ca_coord.append(ca.get_coord())
ca_coord_array = pd.DataFrame(ca_coord)

#-- Create combinations for n X n matrix
import itertools
v_list = []
for i in combinations(ca_coord,2):
    v_list.append(i)
    
#-- Calculate Distances 
joint_list = pd.DataFrame(v_list)
split_df1 = pd.DataFrame(joint_list[0].tolist())
split_df2 = pd.DataFrame(joint_list[1].tolist())
k = pd.concat([split_df1, split_df2], axis=1)
c = k.values.tolist()    
dist_c = []
for i in c:
    dist_c.append(pairwise_atom_distance(i))
    
#-- Network Creation
    #-- Create node labels 1
    #-- n = number of CA atoms in protein
n = count
t2 = []
for y in range(n):
    t2.append(y+1)
    
#-- Create unique combinations
comb_list = []
for i in combinations(t2,2):
    comb_list.append(i)
    
#--  Create source and target lists
convert_tuple_to_list = np.array (comb_list)
source = convert_tuple_to_list[:,0]
target = convert_tuple_to_list[:,1]    
    

#-- Create Network based on cutoff-score from a network edgelist
cutoff = input("Enter distance cutoff: ")
links=[]
a = 0
b = 1
count_a = 0
count_b = 0
for i in dist_c:
    if (i<=int(cutoff)):
        links.append(b)
        count_b = count_b+1
    else:
        links.append(a)
        count_a = count_a+1
        
indices = []
for i in range (len(links)):
       if links[i] == 0:
           indices.append(i)
links_copy = links      
new_links_source = np.delete(source,indices,axis=0)
new_links_target = np.delete(target,indices,axis=0)

edgelist = pd.DataFrame(
    {'source': new_links_source,
     'target': new_links_target,
    })

H = nx.Graph()
H = nx.from_pandas_edgelist(edgelist)
figure(figsize=(20,16))
nx.draw_shell(H, with_labels = True)

#-- Create Distance matrix

comb_list_matrix = [p for p in itertools.product(t2, repeat=2)]
source_matrix =  (np.array(comb_list_matrix)[:,0])
target_matrix =  (np.array(comb_list_matrix)[:,1]) 
zdf_matrix = pd.DataFrame(source_matrix)
zdf_matrix["target"] = target_matrix
zdf_matrix=zdf_matrix.rename(columns={0: "source"})


new_dist_matrix = pd.DataFrame( [p for p in itertools.product(ca_coord, repeat=2)])

split_df1_matrix = pd.DataFrame(new_dist_matrix[0].tolist())
split_df2_matrix = pd.DataFrame(new_dist_matrix[1].tolist())
k_matrix = pd.concat([split_df1_matrix, split_df2_matrix], axis=1)
c_matrix = k_matrix.values.tolist()    
dist_c_matrix = []
for i in c_matrix:
    dist_c_matrix.append(pairwise_atom_distance(i))
 
zdf_matrix["weights"] = dist_c_matrix
df_new_matrix = pd.DataFrame(zdf_matrix)
df_new_matrix = df_new_matrix.pivot(index='source',columns='target',values='weights')
df_new_matrix =df_new_matrix.where(~df_new_matrix.isna(), df_new_matrix.T) #reflects half-matrix across diagonal
df_new_matrix.fillna(0, inplace=True)
nn_matrix = df_new_matrix.to_numpy()

from matplotlib import cm

plt.matshow(df_new_matrix,cmap=cm.pink)
plt.show()

#-- Create Protein Contact Map
bl_w_matrix = df_new_matrix>=int(cutoff)
cl_w_matrix = bl_w_matrix.astype(int)

np.fill_diagonal(cl_w_matrix.values, 1)

plt.matshow(cl_w_matrix,cmap=cm.binary)
 

    
    










