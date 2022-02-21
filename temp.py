import os, re, sys, time
import pandas as pd
import numpy as np
import networkx as nx
from matplotlib.pyplot import figure 
import matplotlib.pyplot as plt
import itertools
from itertools import combinations
from Bio.File import as_handle
from matplotlib import cm

#-- Set Path to PDB Files
os.chdir("C:\\Users\\saySa\\OneDrive\\Desktop\\raw codes PCN")


class MMCIF2Dict(dict):    

    def __init__(self, filename):
        self.quote_chars = ["'", '"']
        self.whitespace_chars = [" ", "\t"]
        with as_handle(filename) as handle:
            loop_flag = False
            key = None
            tokens = self._tokenize(handle)
            try:
                token = next(tokens)
            except StopIteration:
                return  # for Python 3.7 and PEP 479
            self[token[0:5]] = token[5:]
            i = 0
            n = 0
            for token in tokens:
                if token.lower() == "loop_":
                    loop_flag = True
                    keys = []
                    i = 0
                    n = 0
                    continue
                elif loop_flag:               
                    if token.startswith("_") and (n == 0 or i % n == 0):
                        if i > 0:
                            loop_flag = False
                        else:
                            self[token] = []
                            keys.append(token)
                            n += 1
                            continue
                    else:
                        self[keys[i % n]].append(token)
                        i += 1
                        continue
                if key is None:
                    key = token
                else:
                    self[key] = [token]
                    key = None

    def _splitline(self, line):
        in_token = False
    
        quote_open_char = None
        start_i = 0
        for (i, c) in enumerate(line):
            if c in self.whitespace_chars:
                if in_token and not quote_open_char:
                    in_token = False
                    yield line[start_i:i]
            elif c in self.quote_chars:
                if not quote_open_char and not in_token:
                    quote_open_char = c
                    in_token = True
                    start_i = i + 1
                elif c == quote_open_char and (
                    i + 1 == len(line) or line[i + 1] in self.whitespace_chars
                ):
                    quote_open_char = None
                    in_token = False
                    yield line[start_i:i]
            elif c == "#" and not in_token:
                return
            elif not in_token:
                in_token = True
                start_i = i
        if in_token:
            yield line[start_i:]
        if quote_open_char:
            raise ValueError("Line ended with quote open: " + line)

    def _tokenize(self, handle):
        empty = True
        for line in handle:
            empty = False
            if line.startswith("#"):
                continue
            elif line.startswith(";"):

                token_buffer = [line[1:].rstrip()]
                for line in handle:
                    line = line.rstrip()
                    if line.startswith(";"):
                        yield "\n".join(token_buffer)
                        line = line[1:]
                        if line and not line[0] in self.whitespace_chars:
                            raise ValueError("Missing whitespace")
                        break
                    token_buffer.append(line)
                else:
                    raise ValueError("Missing closing semicolon")
            yield from self._splitline(line.strip())
        if empty:
            raise ValueError("Empty file.")
            
#------------------------------------------------------------------------#  
#-- NOTES

'''  Themes:     1. MMCIF to Dict
                 2. Read PDB Dataset from text files
                 2. Extract ATOM information
                 3. PickCA
                 4. Dist_Calc
                 5. Network Creation
                 6. Distance and Adjacency matrix  
                 
 '''  
      
def read_textfile():
    with open('pdb_files.txt','r') as text:
        textfiles = text.readlines()
        suffix_cif = ".cif"
        pdbids=[]
        for element in textfiles:
            pdbids.append(element.strip())
        pdbids = ["{}{}".format(i,suffix_cif) for i in pdbids]
    return pdbids

cif_files = list(read_textfile())

def pairwise_atom_distance(dist_list = [], *args):
    x = list(dist_list)
    x1 = x[0]
    y1 = x[1]
    z1 = x[2]
    x2 = x[3]
    y2 = x[4]
    z2 = x[5]
    import math
    distance = (math.sqrt(((x2-x1)**2) + ((y2-y1)**2) + ((z2-z1)**2)))
    return(distance)

atom_type           = []
atom_site           = []
atom_symbol         = []
atom_fullsymbol     = []
atom_altloc         = []
residue_name        = []
residue_id          = []
insertion_code      = []
x_coord             = []
y_coord             = []
z_coord             = []
occupancy           = []
b_factor            = []
model_number        = [] 
   
file = '4f3l.cif'
atom_type = MMCIF2Dict(file) ['_atom_site.group_PDB']
atom_site = MMCIF2Dict(file) ['_atom_site.id']
atom_symbol = MMCIF2Dict(file) ['_atom_site.type_symbol']
atom_fullsymbol = MMCIF2Dict(file) ['_atom_site.label_atom_id']
atom_altloc = MMCIF2Dict(file) ['_atom_site.label_alt_id']
residue_name = MMCIF2Dict(file) ['_atom_site.label_comp_id']
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
data_cif[13] = data_cif[13].astype(int)
data_cif = data_cif.drop(data_cif[data_cif[13] >1].index)  
    #---   remove altloc 
data_cif = data_cif [~data_cif [4].isin(['B','C','D','E'])]
    
    #---   remove insertion codes
data_cif = data_cif [~data_cif [7].isin(['B','C','D','E'])]

    #--- Pick CA
pick_CA =  data_cif[data_cif[3] == 'CA']
    
    #--- Pick XYZ coordinates
xyz_coord = pick_CA[pick_CA.columns[8:11]]
ca_coord = list(zip(*map(xyz_coord.get, xyz_coord)))
    
    
    #-- Create combinations for n X n matrix

v_list = []
for i in combinations(ca_coord,2):
    v_list.append(i)
        
    #-- Calculate Distances 
joint_list = pd.DataFrame(v_list)
split_df1 = pd.DataFrame(joint_list[0].tolist())
split_df2 = pd.DataFrame(joint_list[1].tolist())
k = pd.concat([split_df1, split_df2], axis=1)
c =  k.values.tolist() 
h = []
for i in c:
    
    h.append(list(map(float, i)))
 

        
dist_c = []
for i in h:
    dist_c.append(pairwise_atom_distance(i))
    
#-- Network Creation
    #-- Create node labels 1
    #-- n = number of CA atoms in protein
n = len(ca_coord)
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

h1 = []
for i in c_matrix:   
    h1.append(list(map(float, i)))
    
dist_c_matrix = []
for i in h1:
    dist_c_matrix.append(pairwise_atom_distance(i))
 
zdf_matrix["weights"] = dist_c_matrix
df_new_matrix = pd.DataFrame(zdf_matrix)
df_new_matrix = df_new_matrix.pivot(index='source',columns='target',values='weights')
df_new_matrix =df_new_matrix.where(~df_new_matrix.isna(), df_new_matrix.T) #reflects half-matrix across diagonal
df_new_matrix.fillna(0, inplace=True)
nn_matrix = df_new_matrix.to_numpy()


labels = pick_CA[6] 


from matplotlib import cm
plt.matshow(df_new_matrix,cmap=cm.pink)
plt.show()

#-- Create Protein Contact Map
bl_w_matrix = df_new_matrix>=int(cutoff)
cl_w_matrix = bl_w_matrix.astype(int)

np.fill_diagonal(cl_w_matrix.values, 1)

plt.matshow(cl_w_matrix,cmap=cm.binary)
   
        
    

