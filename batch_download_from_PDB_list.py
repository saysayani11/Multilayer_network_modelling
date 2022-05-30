#-- DOWNLOAD MMCIF FILES FROM LIST
import os

#-- Set Path
os.chdir('/home/sayantoni/Desktop/Multilayer_Network_Modeling/complex_data')

from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
parser = MMCIFParser()

#-- READ LIST OF PROTEIN COMPLEXES

def read_textfile():
    with open('protein_complex_data_duplicates_removed','r') as text:
        textfiles = text.readlines()
        suffix_cif = ".cif"
        pdbids = []
        mmcif_formatted = []
        for element in textfiles:
            pdbids.append(element.strip())
            
    return pdbids
pdbids  = list(read_textfile())

from pypdb import *
for i in pdbids:
    
    pdb_file = get_pdb_file(i, filetype='cif', compression=False)
    text_file = open(i+".cif", "w")
    text_file.write(pdb_file)
    text_file.close()