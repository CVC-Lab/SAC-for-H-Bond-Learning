# Standard imports
import numpy as np
import os

# Imports to parse MD simulation files and PDB files
import MDAnalysis as mda
from MDAnalysis.analysis.bat import BAT
import nglview as nv

from Bio.PDB import PDBParser, Select, PDBIO

# Dictionary for converting three letter residue abbreviations to 1 letter code 
res_3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
     
# Converts the atomsitic positions to internal coordinates using mdanalysis
# In interanl coordinates, each atom is defined by
# 1. Bond Length (defined by 2 atoms)
# 2. Bond Angle (defined by 3 atoms)
# 3. Dihedral Angle (defined by 4 atoms)
# 
# Returns a N x 4 numpy matrix representing the internal coordinates 
def pdb_to_intcoords(psf, pdb):
    u = mda.Universe(psf, pdb)
    # Select all atoms associated with a protein
    protein_residues = u.select_atoms("protein")
    intern = BAT(protein_residues)
    intern.run()
    return intern

# Returns a NGLView object to visualize a protein
# psf (str): path to the psf file
# coords (str): path to pdb file or traj file (dcd, xtc,...)
def visualize_protein(psf, coords, default=None, default_representation=False):
    u = mda.Universe(psf, coords)
    # Select all atoms associated with a protein
    protein_residues = u.select_atoms("protein")
    w = nv.show_mdanalysis(protein_residues, default=default, default_representation=default_representation)
    return w

# Downloads a given pdb id from the PDB database to output file name
# Default output is the current directory
# pdb_id (str): PDB id
def download_pdb (pdb_id, output_dir=None):
    pdb_id = pdb_id.upper() + ".pdb"
    # Default output is current working directory
    if (output_dir == None):
        output_dir = "./"
    output = os.path.join(output_dir, pdb_id.lower())
    # if the output file does not exists, download the pdb file
    if not os.path.isfile(output):
        os.system("curl -L https://files.rcsb.org/download/{}.gz --output {}.gz".format(pdb_id, output))
        os.system ("gunzip {}.gz".format(output))
    else:
        print ("The file already exists")
    return 

# Grabs coordinates from a file and returns a numpy array of coords.
# Optionally save coords to pdb or npy file
# 
# coord_path (str): path to the trajectory or pdb file
# top_path (str): path to the topology file
# file_type (str): type of coord file (dcd or pdb)
# save_pdbs (bool): if True, saves the coordinates to a pdb file
# Return: N x 3 Numpy array of coords
def get_coords(coord_path, top_path, file_type="dcd", save_pdbs=False, save_np=False, np_file="prot_coords.npy"):
    u = mda.Universe(top_path, coord_path)
    protein = u.select_atoms("protein")
    result = []
    
    # Case if the file is a trajectory file
    if file_type == "dcd":
        # Add the postions of the protein at each timestep to result
        for ts in u.trajectory:
            print(ts)
            result.append(ts.positions)
        # Option to save pdbs of each trajectory
        if (save_pdbs):
            with mda.Writer("protein.pdb", protein.n_atoms) as W:
                for ts in u.trajectory:
                    W.write(protein)
    
    # Case if the file is a single pdb file
    elif file_type == "pdb":
        result.append(u.positions)

    # Return a numpy array of type float64
    result = np.array(result).astype(np.float64)
    if save_np == True:
        np.save(np_file, result)
    return result

# Extract chains from a PDB File and save them as separate PDB files
# 
# file_path (str): file path of the PDB
# chains (array of str): chain ids that we want; default behavior is extract all chains
# output_dir (str): directory to save the extracted chains default behavior is dir of file_path
def extract_chains (file_path, chains=None, output_path=None):
    # Default output path is directory of input file
    if (output_path == None):
        output_path = os.path.dirname(file_path)
    parser = PDBParser (PERMISSIVE=False)
    pdb_id = os.path.splitext(os.path.split(file)[1])[0]
    structure = parser.get_structure (pdb_id, file_path)
    if chains == None:
        chains = [each.id for each in structure.get_chains()]
    # Extract each chain from the pdb file 
    for chain in chains:
        pdb_chain_file = output_path + "/" + pdb_id + "_" + chain.upper() + ".pdb"
        pdb_io = PDBIO()
        pdb_io.set_structure(structure)
        pdb_io.save ('{}'.format(pdb_chain_file), ChainSelect(chain))
    return

# Class used to help select chains
class ChainSelect (Select):
    def __init__(self, chain):
        self.chain = chain
    def accept_chain(self, chain):
        if chain.get_id() == self.chain:
            return 1
        else:          
            return 0

# Receives a pdb file and converts it into a graph representation
# Nodes are atoms and bonds are edges
# Returns a Feature Matrix "X" and an Adjacency Matrix "A"
# TODO STILL HAVE TO IMPLEMENT
def pdb_to_graph(pdb, parm, is_pqr=False):
    # Try to parse the pdb file
    try:
        parser = PDBParser(PERMISSIVE=False, is_pqr=is_pqr)
        structure = parser.get_structure ("id", pdb)
    except:
        print ("There was a parsing error for pdb file at ", pdb)
        return None
    if (is_pqr == False):
        universe = mda.Universe(pdb, parm)
    return X, A

# Converts PDB files to PQR files 
# Missing residues are added with MODELLER
# Missing atoms are added with PDB2PQR 3.1
# Protonation states are assigned with PROPKA3
# Charges and Atom radii are assigned with the AMBER Forcefield FF19SB
# TODO STILL HAVE TO IMPLEMENT
def pdb_to_pqr(pdb):  
    return 

if __name__ == '__main__':
    file = "/home/conrad/Oct-GP/Learning-Viral-Assembly-Pathways-with-RL-/data/6gyp/6gyp.pdb"
    dcd_file = "/mnt/c/Users/conra/CHARMM PDB Sims/1acb_b_rmin_full.dcd"
    pdb_file = "/mnt/c/Users/conra/CHARMM PDB Sims/1acb_b_rmin.pdb"
    top_file = "/mnt/c/Users/conra/CHARMM PDB Sims/1acb_b_rmin.psf"
    coors = get_coords(pdb_file, top_file, save_np=True, np_file="/home/conradli/Learning-Viral-Assembly-Pathways-with-RL-/data/1acb/1acb_coords_pdb.npy")
    print(coors)
    print(coors.shape)
    output_path="/home/conrad/Oct-GP/Learning-Viral-Assembly-Pathways-with-RL-/data/6gyp/chains"
    print("No compile errors")