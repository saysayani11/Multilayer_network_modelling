import os
os.chdir("/mnt/808e843c-fb08-404a-9f45-efd817d81a77/Sayantoni/networks/cif_files_1758")
from pathlib import Path
with open('1758_dataset_19_04_2022','r') as text:
    textfiles = text.readlines()
   
    pdbids=[]
    for element in textfiles:
        pdbids.append(element.strip())
    
os.chdir("/mnt/808e843c-fb08-404a-9f45-efd817d81a77/Sayantoni/networks/results")
count = 0
for items in pdbids:
    filename = str(items)
    directory = "/mnt/808e843c-fb08-404a-9f45-efd817d81a77/Sayantoni/networks/results/" + filename
    os.chdir(directory)
    sub_directory = directory + "/bond_percolation_100/protein_structure_networks/"
    for f in os.listdir(sub_directory):
        os.remove(os.path.join(sub_directory, f))       
    print("removed")
