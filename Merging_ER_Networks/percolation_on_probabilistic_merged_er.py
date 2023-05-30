import os
import numpy as np
import random
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib import cm
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from scipy.spatial.distance import pdist, squareform
from itertools import combinations
from collections import Counter
plt.rcParams.update({'figure.max_open_warning': 0})

import matplotlib.pyplot as plt
import random

def generate_er_network(n, m):
    G = nx.empty_graph(n)
    while not nx.is_connected(G):
        possible_edges = list(nx.non_edges(G))
        random_edge = random.choice(possible_edges)
        G.add_edge(*random_edge)
    while G.number_of_edges() < m:
        possible_edges = list(nx.non_edges(G))
        random_edge = random.choice(possible_edges)
        G.add_edge(*random_edge)
    return G

def select_interactors(G, fraction):
    num_interactors = int(fraction * len(G.nodes))
    return random.sample(G.nodes, num_interactors)

def plot_graph(G, node_colors):
    plt.figure(dpi=300)
    nx.draw(G, node_color=node_colors, with_labels=True)
    plt.show()

n = 24
m = 22
interactor_fraction_G1 = 0.46  # interactor fraction for G1
interactor_fraction_G2 = 0.59  # interactor fraction for G2
interactor_edge_prob = 0.3
G1 = generate_er_network(n, m)
G2 = generate_er_network(n, m)
G2 = nx.relabel_nodes(G2, {i: i + n for i in G2.nodes()})
G1_interactors = select_interactors(G1, interactor_fraction_G1)  # use interactor_fraction_G1 for G1
G2_interactors = select_interactors(G2, interactor_fraction_G2)  # use interactor_fraction_G2 for G2
G3 = nx.union(G1, G2)

for node1 in G1_interactors:
    for node2 in G2_interactors:
        if random.random() < interactor_edge_prob:
            G3.add_edge(node1, node2)




merged_nodes = len(list(G3.nodes))
merged_edges = len(list(G3.edges))


temp = list(G3.edges())
edgelist = pd.DataFrame(temp, columns = ['source', 'target'])
# temp = list(zip(random_sourcelist, random_targetlist))    
# edgelist = pd.DataFrame(temp, columns = ['source', 'target'])
    

indices = edgelist.index.values
H = nx.Graph()
largest_cc = []
for i in range(1000):
    df_len = len(edgelist)      
    remove_frac = int((i/1000)*df_len)
    drop_indices = np.random.choice((indices), remove_frac, replace=False)
    df_subset  = edgelist.drop(drop_indices) 
    H= nx.from_pandas_edgelist(df_subset)     
    largest_cc.append(len(max(nx.connected_components(H), key=len)))
print("done")
raw_y_data = list(largest_cc)

ydata= (raw_y_data - (np.min(raw_y_data))) / (np.max(raw_y_data) - np.min(raw_y_data))
xdata = np.arange(0, 1, 0.001).tolist()

df = pd.DataFrame({'frac_of_removed_links': xdata,
                          'size_of_largest_component_normalized': ydata
                         })
    
