from Bio import AlignIO
import numpy as np
import math
import pandas as pd
from Bio.PDB import PDBParser, PDBIO, Select
import os
os.chdir(r"C:\Users\saySa\OneDrive\Desktop\task_1\MSA\alignment_files\all_alignment")


# Read the MSA
alignment = AlignIO.read("all_18_alignment.fasta", "fasta")


# Extract the sequence for 5zm0
sequence_5zm0 = None
for record in alignment:
    if record.id == '5zm0':  # Adjust the ID if it's different in your alignment
        sequence_5zm0 = record.seq
        break

if not sequence_5zm0:
    raise ValueError("5zm0 sequence not found in the alignment.")

# Calculate Shannon Entropy for each position
def shannon_entropy(column):
    """Calculate Shannon entropy for a given alignment column."""
    unique_base_counts = {char: column.count(char) for char in set(column)}
    total_bases = len(column)
    entropy = sum([-1 * (count/total_bases) * math.log2(count/total_bases) 
                   for count in unique_base_counts.values()])
    return entropy

entropies = [shannon_entropy(column) for column in zip(*alignment)]

# Create the dataframe
df = pd.DataFrame({
    'Sequence_5zm0': sequence_5zm0,
    'Shannon_Entropy': entropies
})

# Filter out rows with gaps in Sequence_5zm0
filtered_df = df[df['Sequence_5zm0'] != '-']

# Save the filtered dataframe to an Excel file
filtered_df.to_excel("filtered_dataframe_of_entropies.xlsx", index=False)


# Extract Shannon entropy values
residue_sequence = list(filtered_df['Sequence_5zm0'])
residue_entropies = list(filtered_df['Shannon_Entropy'])

# Parse the 5zm0 structure from the local file
parser = PDBParser()
structure = parser.get_structure("5zm0", "5zm0.pdb")


# Let's use a direct position-to-entropy map instead of residue names
residue_position_to_entropy = dict(enumerate(residue_entropies, start=1))  # Starts from position 1 for PDB

# Assign B-factors in the PDB structure
for model in structure:
    for chain in model:
        for residue in chain:
            res_id = residue.id[1]  # get the residue position from the PDB structure
            if res_id in residue_position_to_entropy:
                entropy_value = residue_position_to_entropy[res_id]
                for atom in residue:
                    atom.set_bfactor(entropy_value)

# Save the modified structure
io = PDBIO()
io.set_structure(structure)
io.save("5zm0_with_entropy.pdb")


