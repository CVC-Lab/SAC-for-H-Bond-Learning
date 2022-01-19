#%%
# Imports
import os
import sys
# ML Libraries
import matplotlib.pyplot as plt
import torch
import torch.nn as nn
import torch.functional as F
# Protein Analysis Libaries
import MDAnalysis as mda
import nglview as nv
from nglview.datafiles import PDB, XTC
from pyrosetta import *
from pyrosetta.toolbox import *
from pyrosetta.rosetta.core.scoring import ScoreType 
# Custom Libraries
from mol_processors.pdb_processor import download_pdb, get_coords, pdb_to_intcoords, visualize_protein
from mol_processors.RosettaProtein import RosettaProt
%pylab inline
%load_ext autoreload 
%autoreload 2
# %%
# Initialize Rosetta
init()
# Define some constants
PDB_DIR = "./data/1BDD/"
OCTREE_DIR = "./Octree"
# Define some parameters
pdb_id = "1bdd" # The pdb id of the input protein
dcd_file = ""
pdb_init_file = PDB_DIR + "1bdd_pnon.pdb"
pdb_out_file = PDB_DIR + "1cq0_pnon-outo.pdb"
psf_file = PDB_DIR + "1cq0_pnon.psf"
# Set flags
download_protein = True # Flag to download protein
test_prot = RosettaProt(pdb_id=pdb_id)
# %%
'''
print("Initialization")
init()
#path_prefix = "/home/conradli/SAC-for-H-Bond-Learning/Octree/FromBU/oct-example/"
dir_prefix = "/home/conradli/SAC-for-H-Bond-Learning/data"
file_prefix = "ala_dip_charmm"
path_prefix = os.path.join(dir_prefix, "alanine_dipeptide/", file_prefix)
param_prefix = os.path.join(dir_prefix, "params/")
print(path_prefix)
# File Paths for param and pdb files
pdbFile = bytes(path_prefix + ".pdb", encoding="ascii")
pdbFixedFile = bytes(path_prefix + "-fixed.pdb", encoding="ascii") # Currently not used
mol2File = bytes(path_prefix +".mol2", encoding="ascii")
psfFile = bytes(path_prefix + ".psf", encoding="ascii")
outnFile = bytes(path_prefix + "-outn.pdb", encoding="ascii") # Output for nonbonded list
outoFile = bytes(path_prefix + "-outo.pdb", encoding="ascii") # Output for octree
prmFile = bytes(param_prefix + "parm_new.prm", encoding="ascii")
rtfFile = bytes(param_prefix + "pdbamino_new.rtf", encoding="ascii")
aprmFile = bytes(param_prefix + "atoms.0.0.6.prm.ms.3cap+0.5ace.Hr0rec", encoding="ascii")
prot = RosettaProt(pdbFile, psfFile, mol2File, prmFile, rtfFile, aprmFile, outnFile, outoFile)
'''

# %%
#print(test_prot.atom_ids)
print(test_prot.get_cart_coords())
print(len(test_prot.cart_coords))
# %%

# %%
# Memory statprint("Memory Used: ", getrusage(RUSAGE_SELF).ru_maxrss / 100000.0, "MB")
# %%