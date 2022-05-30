#-- CLEANING COMPLEXED PROTEINS CIF FILES 

import os
import pandas as pd
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict


#-- Set Path
os.chdir('/home/sayantoni/Desktop/Multilayer_Network_Modeling/complex_data')


#-- READ LIST OF PROTEIN COMPLEXES, REMOVE NMR STRUCTURES AND FETCH RESOLUTION DETAILS OF XRD STRUCTURES
parser = MMCIFParser()

def read_textfile():
    with open('protein_complex_data_duplicates_removed','r') as text:
        textfiles = text.readlines()
        pdbids = []
        for element in textfiles:
            pdbids.append(element.strip())            
    return pdbids
pdbids  = list(read_textfile())

#-- PARSE PDF FILES IN LOOP, FETCH XRD STRUCTURES ALONG WITH RESOLUTION DETAILS
exptl = []
res = []

for i in pdbids:
    structure = parser.get_structure(i, i+".cif")   
    temp1 = structure.header["structure_method"]
    
    if (temp1 == "X-RAY DIFFRACTION"):
       exptl.append(i)      
       temp2 = structure.header["resolution"]
       res.append(temp2)
       
    
#-- EXPORT/SAVE HEADER DATA OF EACH PDB TO EXCEL SHEET

df = pd.DataFrame(   { "PDBID"  :  exptl,
                        "RESOLUTION"  :  res},
                        columns = ["PDBID", "RESOLUTION"])

df.to_excel("resolution details.xlsx")










    
    

