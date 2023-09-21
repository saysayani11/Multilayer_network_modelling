import os
os.chdir(r"C:\Users\saySa\OneDrive\Desktop\task_1\PSNs\networks")
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib import cm
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from scipy.spatial.distance import pdist, squareform
from itertools import combinations
from collections import Counter
from Bio.PDB import MMCIFParser, Superimposer, PDBIO
from scipy.spatial import distance_matrix



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

def superimpose_common_CA_atoms(ref_structure, mob_structure):
    """Superimposes the mobile structure onto the reference structure using common C-alpha atoms."""
    ref_atoms = extract_CA_atoms(ref_structure)
    mob_atoms = extract_CA_atoms(mob_structure)

    common_atoms = min(len(ref_atoms), len(mob_atoms))
    
    super_imposer = Superimposer()
    super_imposer.set_atoms(ref_atoms[:common_atoms], mob_atoms[:common_atoms])
    super_imposer.apply(mob_structure.get_atoms())

    return super_imposer.rms
 
def save_CA_distance_matrices(ref_atoms, mob_atoms, pair_name):
    """Save the CA-CA distance matrix for reference and mobile atoms."""
    
    # Get coordinates for the atoms
    ref_coords = [atom.get_coord() for atom in ref_atoms]
    mob_coords = [atom.get_coord() for atom in mob_atoms]
    
    # Compute distance matrices
    ref_distance_matrix = distance_matrix(ref_coords, ref_coords)
    mob_distance_matrix = distance_matrix(mob_coords, mob_coords)
    
    # Save the distance matrices
    ref_filename = f"{pair_name}_reference_dist_matrix.csv"
    mob_filename = f"{pair_name}_mobile_dist_matrix.csv"
    
    np.savetxt(ref_filename, ref_distance_matrix, delimiter=",")
    np.savetxt(mob_filename, mob_distance_matrix, delimiter=",")
    
    print(f"Saved distance matrices for {pair_name}.")
    
def save_subtraction_matrix(ref_atoms, mob_atoms, pair_name):
    """Calculate and save the subtraction matrix between the ref and mob structure."""
    
    # Get coordinates for the atoms
    ref_coords = [atom.get_coord() for atom in ref_atoms]
    mob_coords = [atom.get_coord() for atom in mob_atoms]
    
    # Compute distance matrices
    ref_distance_matrix = distance_matrix(ref_coords, ref_coords)
    mob_distance_matrix = distance_matrix(mob_coords, mob_coords)
    
    # Calculate subtraction matrix
    subtraction_matrix = ref_distance_matrix - mob_distance_matrix
    
    # Save the subtraction matrix
    subtraction_filename = f"{pair_name}_subtraction_matrix.csv"
    np.savetxt(subtraction_filename, subtraction_matrix, delimiter=",")
    
    print(f"Saved subtraction matrix for {pair_name}.")
    
    
# Define path and get the list of .cif files
path = r"C:\Users\saySa\OneDrive\Desktop\task_1\PSNs\networks"
os.chdir(path)
cif_files = [os.path.join(path, f) for f in os.listdir() if f.endswith('.cif')]

parser = MMCIFParser()

rmsd_values = {}

for idx, file1 in enumerate(cif_files):
    ref_structure = parser.get_structure("reference", file1)
    for file2 in cif_files[idx+1:]:
        mob_structure = parser.get_structure("mobile", file2)
        
        rmsd = superimpose_common_CA_atoms(ref_structure, mob_structure)
        if rmsd is not None:
            pair_name = f"{os.path.basename(file1).split('.')[0]}_vs_{os.path.basename(file2).split('.')[0]}"
            rmsd_values[pair_name] = rmsd
            
            # Extract common CA atoms again and save distance matrices
            ref_atoms = extract_CA_atoms(ref_structure)
            mob_atoms = extract_CA_atoms(mob_structure)
            common_atoms = min(len(ref_atoms), len(mob_atoms))
            
            save_CA_distance_matrices(ref_atoms[:common_atoms], mob_atoms[:common_atoms], pair_name)
            save_subtraction_matrix(ref_atoms[:common_atoms], mob_atoms[:common_atoms], pair_name)

print(rmsd_values)


