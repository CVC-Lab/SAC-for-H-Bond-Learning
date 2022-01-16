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
from pyrosetta.rosetta.core.chemical import AtomTypeSet
from pyrosetta.rosetta.protocols.relax import FastRelax
# Custom Libraries
from mol_processors.pdb_processor import download_pdb, get_coords, pdb_to_intcoords, visualize_protein
from mol_processors.Protein import Prot
from utils import write_array

# Initialize Rosetta
init()
# Define some constants
PDB_DIR = "./data/1BDD/"
OCTREE_DIR = "./Octree"
OUTPUT_DIR = "./data/1BDD/output"
# Define some parameters
pdb_id = "1bdd" # The pdb id of the input protein
dcd_file = ""
pdb_init_file = PDB_DIR + "1bdd.pdb"
pdb_out_file = PDB_DIR + "1cq0_pdb2pqr_charmm-outn.pdb"
psf_file = PDB_DIR + "1bdd_pdb2pqr_charmm.psf"
# Set flags
download_protein = True # Flag to download protein
test_prot = Prot(pdb_id=pdb_id, pdb_file=pdb_init_file, psf_file=psf_file)

test_prot.sample_bbdep_rotamers()
test_prot.pack_prot()
#print(test_prot.get_sidechain_angles(1))
#print(test_prot.get_sidechain_angles(1))