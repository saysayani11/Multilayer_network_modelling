#-- CLEANING COMPLEXED PROTEINS CIF FILES 
import os

#-- Set Path
os.chdir('/home/sayantoni/Desktop/Multilayer_Network_Modeling/complex_data')

from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
parser = MMCIFParser()

#-- READ LIST OF PROTEIN COMPLEXES AND REMOVE NMR STRUCTURES

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


exptl = []
for i in pdbids:
    structure = parser.get_structure(i, i+".cif")
    temp = structure.header["structure_method"]
    exptl.append(temp)

