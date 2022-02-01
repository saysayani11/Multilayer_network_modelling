import Bio.PDB
from Bio.PDB import PDBList
from Bio.PDB import MMCIFParser

get_pdb = PDBList()
structure_input = input("enter PDB ID: ")
get_pdb.retrieve_pdb_file(structure_input, pdir = '.', file_format = 'mmCif')
parser = MMCIFParser(QUIET = True)
data = parser.get_structure(structure_input, structure_input+".cif")  
 
structure = data[0]
target_atom = structure['A'][110]['CA']
atoms  = list(Bio.PDB.Selection.unfold_entities(structure, 'A'))
ns = Bio.PDB.NeighborSearch(atoms)
close_atoms = list(ns.search(target_atom.coord, 5))


