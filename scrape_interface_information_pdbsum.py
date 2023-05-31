import re
import requests

# =============================================================================
#      Scrape information on interface interactions (A-B only) from PDBSum
# =============================================================================
PDB_ID = "2j4d"
link = "https://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetIface.pl?pdb=" + PDB_ID + "2a1j&chain1=A&chain2=B"
f = requests.get(link)
k = f.text


def extract_numbers(string):
    salt_bridges = re.search(r"Number of salt bridges:\s*(\d+)", string)
    hydrogen_bonds = re.search(r"Number of hydrogen bonds:\s*(\d+)", string)
    non_bonded_contacts = re.search(r"Number of non-bonded contacts:\s*(\d+)", string)
    disulphide_bonds = re.search(r"Number of disulphide bonds:\s*(\d+)", string)

    if salt_bridges:
        salt_bridges = int(salt_bridges.group(1))
    else:
        salt_bridges = 0

    if hydrogen_bonds:
        hydrogen_bonds = int(hydrogen_bonds.group(1))
    else:
        hydrogen_bonds = 0

    if non_bonded_contacts:
        non_bonded_contacts = int(non_bonded_contacts.group(1))
    else:
        non_bonded_contacts = 0

    if disulphide_bonds:
        disulphide_bonds = int(disulphide_bonds.group(1))
    else:
        disulphide_bonds = 0

    total = salt_bridges + hydrogen_bonds + non_bonded_contacts + disulphide_bonds
    return total


string = k
print(extract_numbers(string))
