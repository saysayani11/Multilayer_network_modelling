import networkx as nx
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

n = 10
m = 20
interactor_fraction = 0.2
G1 = generate_er_network(n, m)
G2 = generate_er_network(n, m)
G2 = nx.relabel_nodes(G2, {i: i + n for i in G2.nodes()})
G1_interactors = select_interactors(G1, interactor_fraction)
G2_interactors = select_interactors(G2, interactor_fraction)

node_colors_G1 = ["lightblue" if node in G1_interactors else "blue" for node in G1.nodes()]
node_colors_G2 = ["pink" if node in G2_interactors else "red" for node in G2.nodes()]

plot_graph(G1, node_colors_G1)
plot_graph(G2, node_colors_G2)
G3 = nx.union(G1, G2)

for node1 in G1_interactors:
    for node2 in G2_interactors:
        G3.add_edge(node1, node2)

node_colors = ["lightblue" if node in G1_interactors else "blue" if node in G1.nodes() else "pink" if node in G2_interactors else "red" for node in G3.nodes()]
plot_graph(G3, node_colors)
