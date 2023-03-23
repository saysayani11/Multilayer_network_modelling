import networkx as nx
import random
import matplotlib.pyplot as plt

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

# Simulate bond percolation and plot results
num_repeats = 100  # Number of times to repeat simulation at each edge removal probability
frac_range = [i / 100 for i in range(101)]  # Range of fraction of edges removed to simulate
largest_components = []  # List to store size of largest component for each fraction of edges removed

for frac in frac_range:
    largest_component_size = 0
    for i in range(num_repeats):
        # Copy merged network and remove edges with given probability
        G = merged.copy()
        for (u, v) in G.edges():
            if random.random() < frac:
                G.remove_edge(u, v)

        # Calculate size of largest connected component
        largest_cc = max(nx.connected_components(G), key=len)
        largest_cc_size = len(largest_cc)

        # Update largest component size if larger than previous repeats
        if largest_cc_size > largest_component_size:
            largest_component_size = largest_cc_size

    largest_components.append(largest_component_size/n)

# Plot results
plt.plot(frac_range, largest_components, linestyle='dashed', linewidth=0.2, marker='o', markersize=3, color='blue')
plt.xlabel('Fraction of edges removed')
plt.ylabel('Fraction of nodes in largest component')
plt.title('Bond percolation on Merged ER network')
plt.ylim([0, 1])
# plt.figsize=(8,6)
plt.show()
