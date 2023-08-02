import os
import pandas as pd
from Bio import SeqIO
os.chdir(r"C:\Users\saySa\OneDrive\Desktop\task_1\sequences_19")


#-- Ranges for Antenna Binding Domain (ABD) and FAD Binding Domain (FAD-BD)
#-- Data from InterProScan Database v5
ABD_ranges = [ (5, 223), (5, 227), (38, 207), (19, 180), (31, 161), (6, 182), (5, 172), (6, 168), (3, 166), (3, 162), (21, 175), (6, 168), (179, 353), (84, 250), (3, 164), (2, 168), (6, 165), (12, 169), (5, 160) ]
FADBD_ranges = [ (226, 491), (229, 488), (234, 488), (228, 461), (252, 436), (213, 514), (290, 487), (297, 494), (286, 483), (288, 486), (306, 504), (288, 478), (393, 508), (368, 526), (241, 395), (272, 466), (274, 473), (290, 488), (287, 485) ]


def extract_domains(fasta_file):
    ABD_sequences = []
    FADBD_sequences = []

    #-- Read the FASTA file
    for i, record in enumerate(SeqIO.parse(fasta_file, "fasta")):
        ABD_start, ABD_end = ABD_ranges[i]
        FADBD_start, FADBD_end = FADBD_ranges[i]
        
        #-- Extract & append the sequences
        ABD_sequences.append(str(record.seq[ABD_start-1:ABD_end]))
        FADBD_sequences.append(str(record.seq[FADBD_start-1:FADBD_end]))

    #-- Store
    ABD_df = pd.DataFrame({'SL': range(1, 20), 'ABD_Sequences': ABD_sequences})
    FADBD_df = pd.DataFrame({'SL': range(1, 20), 'FADBD_Sequences': FADBD_sequences})

    return ABD_df, FADBD_df

fasta_file = 'sequences_19.fasta'
ABD_df, FADBD_df = extract_domains(fasta_file)
print(ABD_df)
print(FADBD_df)
