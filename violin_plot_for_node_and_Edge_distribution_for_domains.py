import seaborn as sns
import matplotlib.pyplot as plt
import os
os.chdir(r"C:\Users\saySa\OneDrive\Desktop\task_1\PSNs\networks\Degree")


# Sample data for nodes and edges in Antenna Binding Domain and FAD Binding Domain
nodes_antenna = [219, 223, 170, 162, 131, 177, 168, 163, 164, 160, 155, 163, 175, 167, 162, 167, 160, 158]
edges_antenna = [108, 105, 47, 44, 34, 52, 44, 39, 42, 43, 46, 39, 50, 64, 44, 45, 47, 35]

nodes_fad = [266, 260, 239, 234, 185, 302, 198, 198, 199, 189, 191, 104, 48, 155, 199, 232, 199, 103]
edges_fad = [115, 113, 109, 99, 58, 124, 74, 56, 75, 76, 78, 70, 42, 68, 49, 103, 66]

# Set DPI to 300
plt.figure(dpi=300)

# Use the "magma" color palette
colors = sns.color_palette("viridis", n_colors=2)

# Create a figure with two subplots
fig, axes = plt.subplots(1, 2, figsize=(12, 6))

# Violin plot for Nodes with the "magma" color palette
sns.violinplot(data=[nodes_antenna, nodes_fad], ax=axes[0], inner="stick", palette=colors)
axes[0].set_xticklabels(['Antenna Binding Domain', 'FAD Binding Domain'])
axes[0].set_ylabel('Number of Nodes')
axes[0].set_title('Distribution of Nodes')

# Increase axes ticks
axes[0].tick_params(axis='both', which='both', length=6, width=2, labelsize=10)

# Violin plot for Edges with the "magma" color palette
sns.violinplot(data=[edges_antenna, edges_fad], ax=axes[1], inner="stick", palette=colors)
axes[1].set_xticklabels(['Antenna Binding Domain', 'FAD Binding Domain'])
axes[1].set_ylabel('Number of Edges')
axes[1].set_title('Distribution of Edges')

# Increase axes ticks
axes[1].tick_params(axis='both', which='both', length=12, width=2, labelsize=10)

plt.tight_layout()



# Save the plots with DPI 300
plt.savefig('violin_plots.png', dpi=300)
plt.show()
