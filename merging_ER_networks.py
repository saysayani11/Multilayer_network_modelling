import networkx as nx
import random

# Generate two ER networks using G(n,m) model
n = 50  # Number of nodes
m = 2  # Number of edges for each node
er1 = nx.gnm_random_graph(n, m * n // 2)
er2 = nx.gnm_random_graph(n, m * n // 2)

# Randomly select some nodes as interactors
num_interactors = 10
interactors1 = random.sample(er1.nodes(), num_interactors)
interactors2 = [n + i for i in interactors1]  # To get the same nodes in er2

# Merge the two networks
merged = nx.Graph()
merged.add_nodes_from(er1.nodes())
merged.add_nodes_from(er2.nodes())
merged.add_edges_from(er1.edges())
merged.add_edges_from(er2.edges())
for i, j in zip(interactors1, interactors2):
    merged.add_edge(i, j)
    non_interactors = set(er1.nodes()).difference(set(interactors1))
    non_interactor = random.sample(non_interactors, 1)[0]
    merged.add_edge(i, non_interactor)

# Print the initial and final node and edge counts
print("ER1: nodes = {}, edges = {}".format(er1.number_of_nodes(), er1.number_of_edges()))
print("ER2: nodes = {}, edges = {}".format(er2.number_of_nodes(), er2.number_of_edges()))
print("Merged: nodes = {}, edges = {}".format(merged.number_of_nodes(), merged.number_of_edges()))
