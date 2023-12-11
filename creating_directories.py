import os
os.chdir("/mnt/808e843c-fb08-404a-9f45-efd817d81a77/Sayantoni/networks/cif_files_1758")
  
with open('1758_dataset_19_04_2022','r') as text:
    textfiles = text.readlines()
   
    pdbids=[]
    for element in textfiles:
        pdbids.append(element.strip())
    
os.chdir("/mnt/808e843c-fb08-404a-9f45-efd817d81a77/Sayantoni/networks/results")

for items in pdbids:    
    os.mkdir(items)
