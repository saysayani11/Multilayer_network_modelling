from Bio.PDB import PDBParser, MMCIFIO
import os
os.chdir(r"C:\Users\saySa\OneDrive\Desktop\task_1\PSNs")

pdb_parser = PDBParser(QUIET=True)
mmcif_io = MMCIFIO()

with open('str_files.txt', 'r') as protein_file:
    for line in protein_file:
        protein_name = line.strip()
        

        pdb_path = f"./{protein_name}.pdb"
        mmcif_path = f"./{protein_name}.cif"

        structure = pdb_parser.get_structure(protein_name, pdb_path)
        mmcif_io.set_structure(structure)
        mmcif_io.save(mmcif_path)
