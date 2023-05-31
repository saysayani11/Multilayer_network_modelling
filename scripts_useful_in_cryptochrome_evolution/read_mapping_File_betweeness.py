import networkx as nx
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser, is_aa
import os
import pandas as pd
#--- Set Path
os.chdir(r"C:\Users\saySa\OneDrive\Desktop")

# Read the mapping file and store the betweenness centrality values
mapping_file = 'mapping.xlsx'
df = pd.read_excel(mapping_file)
df.dropna(inplace=True)
betweenness_values = dict(zip(df['Residue'].astype(str), df['Betweenness'].astype(float)))

# Load the PDB structure using Biopython
pdb_id = '3umv'
pdb_file = f'{pdb_id}.pdb'

parser = PDBParser()
structure = parser.get_structure(pdb_id, pdb_file)

# Create a NetworkX graph
graph = nx.Graph()