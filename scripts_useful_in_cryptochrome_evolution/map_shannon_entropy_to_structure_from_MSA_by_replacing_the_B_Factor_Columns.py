from Bio import AlignIO
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select
import os
os.chdir(r"C:\Users\saySa\OneDrive\Desktop\task_1\MSA\alignment_files\all_alignment")

# Calculate the Shannon entropy of a column in the MSA
def shannon_entropy(col):
    unique_base, counts = np.unique([base for base in col if base != '-'], return_counts=True)
    prob = counts / float(sum(counts))
    entropy = -sum(p * np.log2(p) for p in prob)
    return entropy

# Read the MSA
alignment = AlignIO.read("all_18_alignment.fasta", "fasta")
msa_cols = np.array([list(rec) for rec in alignment], np.character)

# Calculate Shannon entropy for each column in the MSA
entropies = [shannon_entropy(msa_cols[:, i]) for i in range(msa_cols.shape[1])]

# Load the 5zm0 structure
parser = PDBParser()
structure = parser.get_structure("5zm0", "5zm0.pdb")

# Map Shannon entropy to B-factor
for model in structure:
    for chain in model:
        for i, residue in enumerate(chain):
            for atom in residue:
                # Use the Shannon entropy as the B-factor
                if i < len(entropies):
                    atom.set_bfactor(entropies[i])

# Save the modified structure
io = PDBIO()
io.set_structure(structure)
io.save("5zm0_entropy.pdb")
