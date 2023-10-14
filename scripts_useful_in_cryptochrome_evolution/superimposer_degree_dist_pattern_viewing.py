import os
import numpy as np
import pandas as pd
from scipy.stats import entropy
from scipy.spatial import distance
from Bio.PDB import Superimposer
from Bio.PDB import MMCIFParser
from itertools import combinations
from Bio.PDB import MMCIFParser, Superimposer, PDBIO
from scipy.spatial.distance import pdist, squareform
from collections import Counter
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.mplot3d import Axes3D
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
plt.rcParams.update({'figure.max_open_warning': 0, 'figure.dpi': 300})
from Bio.PDB import MMCIFParser
parser = MMCIFParser(QUIET=True)

os.chdir(r"C:\Users\saySa\OneDrive\Desktop\task_1\PSNs\networks")

def read_textfile():
    with open('test.txt','r') as text:
        textfiles = text.readlines()
        suffix_cif = ".cif"
        pdbids=[]
        for element in textfiles:
            pdbids.append(element.strip())
        pdbids = ["{}{}".format(i,suffix_cif) for i in pdbids]
    return pdbids
cif_files = list(read_textfile())

def fetch(file):
    mmcif_dict = MMCIF2Dict(file)
    data_cif = pd.DataFrame({
        'atom_type': mmcif_dict['_atom_site.group_PDB'],
        'atom_site': mmcif_dict['_atom_site.id'],
        'atom_symbol': mmcif_dict['_atom_site.type_symbol'],
        'atom_fullsymbol': mmcif_dict['_atom_site.label_atom_id'],
        'atom_altloc': mmcif_dict['_atom_site.label_alt_id'],
        'residue_name': mmcif_dict['_atom_site.label_comp_id'],
        'atom_assym': mmcif_dict['_atom_site.label_asym_id'],
        'residue_id': mmcif_dict['_atom_site.label_seq_id'],
        'insertion_code': mmcif_dict['_atom_site.pdbx_PDB_ins_code'],
        'x_coord': mmcif_dict['_atom_site.Cartn_x'],
        'y_coord': mmcif_dict['_atom_site.Cartn_y'],
        'z_coord': mmcif_dict['_atom_site.Cartn_z'],
        'occupancy': mmcif_dict['_atom_site.occupancy'],
        'b_factor': mmcif_dict['_atom_site.B_iso_or_equiv'],
        'model_number': mmcif_dict['_atom_site.pdbx_PDB_model_num']
    })
    
    # Filter data
    data_cif = data_cif[
            ~data_cif['atom_type'].isin(['HETATM', 'TER']) &
            (data_cif['model_number'].astype(int) == 1) &
            ~data_cif['atom_altloc'].isin(['B', 'C', 'D', 'E']) &
            ~data_cif['insertion_code'].isin(['B', 'C', 'D', 'E']) &
            (data_cif['atom_assym'] == 'A') &
            (data_cif['atom_fullsymbol'] == 'CA')
    ]
    
    data_cif = data_cif[['residue_id', 'x_coord', 'y_coord', 'z_coord', 'atom_site', 'atom_assym', 'residue_name']]
    data_cif.columns = ['residue_no', 'x_coord', 'y_coord', 'z_coord', 'ca_atom_number', 'chain', 'residue_name']
    
    return data_cif

def extract_CA_atoms(structure):
    """Extracts C-alpha atoms from the structure."""
    atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                try:
                    atoms.append(residue["CA"])
                except:
                    pass
    return atoms


def _superimpose(file):
    """Superimpose the mobile structure onto the reference structure using common C-alpha atoms and save the superimposed structure."""
    
    # Parse the structures
    ref_structure = MMCIFParser().get_structure('reference', "3zxs.cif")
    mob_structure = MMCIFParser().get_structure('mobile', file)

    # Extract C-alpha atoms
    ref_atoms = extract_CA_atoms(ref_structure)
    mob_atoms = extract_CA_atoms(mob_structure)
    common_atoms = min(len(ref_atoms), len(mob_atoms))

    # Perform the superimposition
    super_imposer = Superimposer()
    super_imposer.set_atoms(ref_atoms[:common_atoms], mob_atoms[:common_atoms])
    super_imposer.apply(mob_structure.get_atoms())
    
    # Save the superimposed structure as a PDB file
    pdb_io = PDBIO()
    pdb_io.set_structure(mob_structure)
    pdb_io.save(f"superimposed_{os.path.splitext(file)[0]}.pdb")
    
    return mob_structure, super_imposer.rms

def _plot(file):
    """Plot the superimposed reference and mobile structures using matplotlib 3D with edges and specific aesthetics."""

    # Superimpose the structures
    mob_structure, rms = _superimpose(file)
    
    # Extract C-alpha atoms' coordinates from the superimposed mobile structure and reference structure
    ref_atoms = extract_CA_atoms(MMCIFParser().get_structure('reference', "3zxs.cif"))
    mob_atoms = extract_CA_atoms(mob_structure)
    
    ref_coords = [atom.get_coord() for atom in ref_atoms]
    mob_coords = [atom.get_coord() for atom in mob_atoms]

    # Create edges based on the 8 Ångström cutoff
    ref_edges = [(i, j) for i, coord1 in enumerate(ref_coords) for j, coord2 in enumerate(ref_coords) 
                 if i != j and distance.euclidean(coord1, coord2) < 8]
    mob_edges = [(i, j) for i, coord1 in enumerate(mob_coords) for j, coord2 in enumerate(mob_coords) 
                 if i != j and distance.euclidean(coord1, coord2) < 8]

    # Prepare the data for plotting
    ref_x, ref_y, ref_z = zip(*ref_coords)
    mob_x, mob_y, mob_z = zip(*mob_coords)

    # Create a new figure
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot the atoms for both structures with reduced marker size
    ax.scatter(ref_x, ref_y, ref_z, c='blue', marker='o', s=4, label="3zxs")
    ax.scatter(mob_x, mob_y, mob_z, c='red', marker='o', s=4, label=f"{os.path.splitext(file)[0]}")

    # Plot the edges for both structures 
    for (start, end) in ref_edges:
        ax.plot([ref_coords[start][0], ref_coords[end][0]],
                [ref_coords[start][1], ref_coords[end][1]],
                [ref_coords[start][2], ref_coords[end][2]], 'b-', linewidth=0.3)
        
    for (start, end) in mob_edges:
        ax.plot([mob_coords[start][0], mob_coords[end][0]],
                [mob_coords[start][1], mob_coords[end][1]],
                [mob_coords[start][2], mob_coords[end][2]], 'r-', linewidth=0.3)

    # Remove axes and grid
    ax.axis("off")
    ax.grid(False)

    # Set title and legend
    ax.legend(loc='upper left', frameon=False)
    
    # Display the plot
    plt.show()

def _degree(file):
    """Calculate the degree of all the CA atoms (nodes) in the given structure file."""
    
    # Parse the structure
    structure = MMCIFParser().get_structure('protein', file)

    # Extract C-alpha atoms
    atoms = extract_CA_atoms(structure)
    
    # Extract coordinates for these atoms
    coords = [atom.get_coord() for atom in atoms]
    
    # Calculate pairwise distance matrix
    dist_matrix = squareform(pdist(coords))
    
    # Find atom pairs (edges) within the 8 Ångströms cutoff
    edges = np.where((dist_matrix > 0) & (dist_matrix < 8))
    
    # Calculate degree for each atom
    degrees = Counter(edges[0])
    
    return degrees

def superimpose_common_CA_atoms(ref_structure, mob_structure):
    """Superimposes the mobile structure onto the reference structure using common C-alpha atoms."""
    ref_atoms = extract_CA_atoms(ref_structure)
    mob_atoms = extract_CA_atoms(mob_structure)

    common_atoms = min(len(ref_atoms), len(mob_atoms))
    
    super_imposer = Superimposer()
    super_imposer.set_atoms(ref_atoms[:common_atoms], mob_atoms[:common_atoms])
    super_imposer.apply(mob_structure.get_atoms())

    return super_imposer.rms

def _colour_by_degree(file):
    # Load the structures
    ref_structure = parser.get_structure("reference", "3zxs.cif")
    mob_structure = parser.get_structure("mobile", file)

    # Superimpose structures
    rms = superimpose_common_CA_atoms(ref_structure, mob_structure)

    ref_atoms = extract_CA_atoms(ref_structure)
    mob_atoms = extract_CA_atoms(mob_structure)

    # Calculate adjacency matrices for both structures
    ref_coords = np.array([atom.get_coord() for atom in ref_atoms])
    mob_coords = np.array([atom.get_coord() for atom in mob_atoms])

    ref_distances = pdist(ref_coords)
    mob_distances = pdist(mob_coords)

    ref_adj_matrix = np.where(squareform(ref_distances) < 8.0, 1, 0)
    mob_adj_matrix = np.where(squareform(mob_distances) < 8.0, 1, 0)

    # Calculate degree for each atom (node)
    ref_degrees_list = np.sum(ref_adj_matrix, axis=1)
    mob_degrees_list = np.sum(mob_adj_matrix, axis=1)

    # Start plotting
    fig = plt.figure(figsize=(12, 6))
    cmap = plt.cm.viridis

    # Reference structure plot
    ax1 = fig.add_subplot(121, projection='3d')
    ref_scatter = ax1.scatter(ref_coords[:,0], ref_coords[:,1], ref_coords[:,2], c=ref_degrees_list, cmap=cmap, s=6)
    ax1.set_title("Reference 3zxs")

    # Mobile structure plot
    ax2 = fig.add_subplot(122, projection='3d')
    mob_scatter = ax2.scatter(mob_coords[:,0], mob_coords[:,1], mob_coords[:,2], c=mob_degrees_list, cmap=cmap, s=6)
    ax2.set_title(f"Superimposed {file.split('.')[0]} (RMSD: {rms:.2f} Å)")

    # Draw edges
    for i, atom_coord in enumerate(ref_coords):
        for j, other_coord in enumerate(ref_coords):
            if ref_adj_matrix[i, j] == 1:
                ax1.plot3D(*zip(*[atom_coord, other_coord]), color=cmap(ref_degrees_list[i]/max(ref_degrees_list)), linewidth=0.3)

    for i, atom_coord in enumerate(mob_coords):
        for j, other_coord in enumerate(mob_coords):
            if mob_adj_matrix[i, j] == 1:
                ax2.plot3D(*zip(*[atom_coord, other_coord]), color=cmap(mob_degrees_list[i]/max(mob_degrees_list)), linewidth=0.3)

    # Adjusting the visualization
    for ax in [ax1, ax2]:
        ax.axis("off")
        ax.grid(None)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])

    # Colorbars
    cb1 = fig.colorbar(ref_scatter, ax=ax1, orientation='vertical', fraction=0.03, pad=0.1)
    cb1.set_label('Degree')
    cb2 = fig.colorbar(mob_scatter, ax=ax2, orientation='vertical', fraction=0.03, pad=0.1)
    cb2.set_label('Degree')

    plt.tight_layout()
    plt.show()
    
    
for i in cif_files:
    data = fetch(i)
    rmsd = _superimpose(i)
    print(f"RMSD for {i}: {rmsd}")
    # _plot(i)
    _degree(i)
    _colour_by_degree(i)



    
    
