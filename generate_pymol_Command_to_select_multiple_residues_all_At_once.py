
#--- Provided list
residues = [
    ('VAL', 7),
    ('GLY', 9),
    ('GLN', 11),
    ('LEU', 12),
    ('GLU', 37),
    ('GLY', 38),
    ('HIS', 43),
    ('LYS', 47),
    ('PHE', 58),
    ('VAL', 100),
    ('THR', 102),
    ('ARG', 103),
    ('PRO', 104),
    ('GLY', 105),
    ('ASP', 106),
    ('PHE', 128),
    ('PHE', 215),
    ('MET', 252),
    ('SER', 266),
    ('ASN', 269),
    ('LEU', 270),
    ('LEU', 273),
    ('PRO', 275),
    ('LEU', 304),
    ('ARG', 311),
    ('TYR', 337),
    ('ALA', 350),
    ('GLN', 353),
    ('THR', 354),
    ('TYR', 359),
    ('ALA', 360),
    ('PHE', 372),
    ('VAL', 383),
    ('TRP', 386),
    ('TYR', 387),
    ('LEU', 388),
    ('VAL', 390),
    ('TYR', 391),
    ('ASP', 393),
    ('ALA', 394),
    ('LEU', 395),
    ('GLU', 396),
    ('TRP', 397),
    ('VAL', 398),
    ('GLU', 399),
    ('ALA', 400),
    ('PRO', 401),
    ('THR', 403),
    ('MET', 406),
    ('SER', 407),
    ('LEU', 415),
]

#--- Initialize command string with the first selection
command_string = f"select my_selection, (resn {residues[0][0]} and resi {residues[0][1]})"

#--- Add the rest of the residues to the selection
for resn, resi in residues[1:]:
    command_string += f" or (resn {resn} and resi {resi})"

print(command_string)
