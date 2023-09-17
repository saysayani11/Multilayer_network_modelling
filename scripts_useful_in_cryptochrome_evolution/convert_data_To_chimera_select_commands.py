# Original data
original_data = "446_TRP 238_ILE 458_CYS 406_HIS 337_TRP 436_ASN 432_PRO 435_LYS 352_HIS 366_GLY 361_CYS 465_PRO 434_LEU 414_VAL 424_GLY 470_ASN 456_ALA 240_ASN 356_ARG 321_ASP 455_GLN 342_MET 332_GLN 489_GLN 319_TYR 428_ARG 430_TYR 351_ILE 453_GLN 334_GLY 320_GLU 450_GLU 468_MET 425_GLN 474_ALA 476_ASP 433_VAL 469_VAL 484_ARG 437_PHE 427_ILE 416_PHE 478_ASN 363_LEU 243_LEU 488_GLU 438_PRO 354_LEU 325_LEU 426_TYR 475_SER 358_ALA 409_ARG 241_SER 485_VAL 449_SER 408_THR 417_GLY 472_LYS 246_THR 452_GLU 336_PRO 429_LYS 462_ARG 324_ARG 482_MET 359_VAL 407_TYR 421_ASP 451_GLU 365_ARG 245_SER 466_PHE 461_GLY 360_ALA 471_HIS 412_CYS 331_ALA 357_HIS 445_PRO 483_ARG 343_THR 345_LEU 244_PRO 467_PRO 440_LYS 411_PHE 340_ALA 239_PRO 344_GLN 444_GLU 454_ARG 480_GLN 248_GLY 353_HIS 487_GLU 459_ILE 431_LEU 370_ILE 464_TYR 318_TRP 367_ASP 327_LYS 362_PHE 328_TRP 419_ARG 339_ASP 341_ILE 460_ILE 323_GLU 413_PRO 473_GLU 422_PRO 330_THR 457_GLY 477_ARG 355_ALA 418_LYS 237_THR 442_ILE 415_ARG 448_ALA 346_ARG 420_THR 463_ASP 242_LEU 447_THR 441_TYR 329_LYS 338_ILE 335_PHE 481_LEU 368_LEU 369_TRP 333_THR 439_THR 486_ARG 410_ILE 423_GLU 326_HIS 364_THR 443_TYR 479_LEU 322_ALA"
# Split the original data into individual entries
entries = original_data.strip().split()

# Reformat the entries into the desired format
data_community_1 = [entry.split("_")[0] for entry in entries]

# Print the reformatted data
print(data_community_1)



# Extract the residue numbers from the data and convert them to a string with commas
residue_numbers = [residue.split("_")[0] for residue in data_community_1]
residue_numbers_str = ', '.join(residue_numbers)  # Add a comma and space between residue numbers

# Create a selection command using the extracted residue numbers
selection_command = f"select :{residue_numbers_str}"

# Store the selection command in a variable
stored_selection = selection_command

# Print the stored selection command (optional)
print("Stored Selection Command:")
print(stored_selection)
