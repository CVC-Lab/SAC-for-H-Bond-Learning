# Constants for processing Dunbrack ramachandran files (these constants are 
# specfic for files found in data/torsion_libs/ramachandran)
RAMA_LINES_PER_RESIDUE = 108864
RAMA_LINES_PER_NEIGHBOR = 5184
RAMA_CHARS_PER_LINE = 65
RAMA_RES_ORDER = {"ALA": 0, "ARG": 1, "ASN": 2, "ASP": 3, "CPR": 4, "CYS": 5,
                 "GLN": 6, "GLU": 7, "GLY": 8, "HIS": 9, "ILE": 10, 
                 "LEU": 11, "LYS": 12,"MET": 13, "PHE": 14, "PRO": 15, 
                 "SER": 16, "THR": 17, "TRP": 18, "TYR": 19, "VAL": 20}
RAMA_NEIGHBOR_ORDER = {"ALA": 0, "ARG": 1, "ASN": 2, "ASP": 3, "CYS": 4,
                 "GLN": 5, "GLU": 6, "GLY": 7, "HIS": 8, "ILE": 9, 
                 "LEU": 10, "LYS": 11,"MET": 12, "PHE": 13, "PRO": 14, 
                 "SER": 15, "THR": 16, "TRP": 17, "TYR": 18, "VAL": 19}