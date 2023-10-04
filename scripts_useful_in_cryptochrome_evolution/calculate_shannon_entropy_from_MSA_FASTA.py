import matplotlib.pyplot as plt
import seaborn as sns
from Bio import AlignIO
import math
import os
os.chdir(r"C:\Users\saySa\OneDrive\Desktop\task_1\MSA\alignment_files\all_alignment")
import pandas as pd

def calculate_entropy(position):
    freq = {}
    for amino_acid in position:
        if amino_acid not in freq:
            freq[amino_acid] = 1
        else:
            freq[amino_acid] += 1
    
    entropy = 0
    for key in freq:
        probability = freq[key] / len(position)
        entropy -= probability * math.log2(probability)
        
    return entropy

def msa_entropy(msa_file, file_format="fasta"):
    alignment = AlignIO.read(msa_file, file_format)
    length = alignment.get_alignment_length()
    entropies = []

    for i in range(length):
        position = [record.seq[i] for record in alignment]
        entropy = calculate_entropy(position)
        entropies.append(entropy)
    
    return entropies

# Usage:
msa_file = "all_18_alignment.fasta"  # Make sure to provide the right filename here
entropies = msa_entropy(msa_file)

# Visualization
plt.figure(figsize=(15, 5))
sns.lineplot(x=range(1, len(entropies)+1), y=entropies, lw=2)
plt.fill_between(range(1, len(entropies)+1), entropies, color='skyblue', alpha=0.4)
plt.title("Variability Across Sequence Positions", fontsize=20)
plt.xlabel("Position", fontsize=16)
plt.ylabel("Shannon Entropy", fontsize=16)
plt.grid(True, which="both", ls="--", c='0.7')

# Highlight regions with high entropy (You can set a custom threshold if needed)
threshold = max(entropies) * 0.8
plt.axhline(y=threshold, color='r', linestyle='--')
plt.text(len(entropies) * 0.8, threshold + 0.05, 'High Variability Threshold', color = 'r')

plt.tight_layout()
plt.show()


df = pd.DataFrame({'Position': range(1, len(entropies)+1), 'Entropy': entropies})
df.to_excel("entropy_values.xlsx", index=False)

