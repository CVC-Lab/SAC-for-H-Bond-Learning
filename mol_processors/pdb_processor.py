# Standard imports
import numpy as np
import os
from Bio.PDB import PDBParser, Select, PDBIO
import nglview as nv

# Dictionary for converting three letter residue abbreviations to 1 letter code 
res_3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

# TODO: implement getting atom indices for a residue and `set_coords`
# TODO: Verify that the PDBParser actually works

class PDBParser:
    """
    Uses biopython to extract protein information from PDBs or PQR files. Currently, there
    can only be one PDBParser per chain.
    """
    def __init__(self, pdb_path, chain_id=None):
        parser = PDBParser()
        self.structure = parser.get_structure("mol", pdb_path)

        # Grab first chain of first model as the protein if not specfied
        if chain_id is None:
            self.chain = self.structure[0]["A"]
        self.total_res = len(self.chain)
        self.residues = list(self.chain.get_residues())
        return
    
    def get_res_atom_indices(self, res_id):
        """
        Grabs the indices of the atoms in a residue
        """
        return

    def get_num_residues(self):
        return self.total_res

    def get_coords(self):
        """
        Extracts the atomic coordinates from the PDB
        """
        coords = []
        for atom in self.chain.get_atoms():
            coords.append(atom.get_coord())
        return np.array(self.coords)
    
    def set_coords(self, new_coords):
        """
        Sets the coordinates of the protein in the PDB Parser. This is mainly to allow user
        to visualize what a protein looks like after a transformation with `visualize_protein`
        """
        return

    def get_charges(self):
        """
        Extracts the charges for each atom. This will only return relevant information if
        the PDB is a PQR file
        """
        charges = []
        for atom in self.chain.get_atoms():
            charges.append(atom.get_charge())
        return np.array(charges)

    def get_radii(self):
        """
        Extracts the atomic radii for each atom. This will only return relevant information if
        the PDB is a PQR file
        """
        radii = []
        for atom in self.chain.get_atoms():
            radii.append(atom.get_radius())
        return np.array(radii)
    
    def visualize_protein(self, default=None, default_representation=False):
        """
        Visualize a protein structure using Nglview
        Args:
            default (str): the default color for the protein
            default_representation (str): default representation for the protein
        
        Returns:
            viz: the NGLWidget object
        """
        viz = nv.show_biopython(self.structure, default=default, default_representation=default_representation)
        return viz

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
