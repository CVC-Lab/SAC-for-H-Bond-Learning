import MDAnalysis as mda
import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.internal_coords import IC_Chain, IC_Residue
from sklearn.metrics import pairwise_distances_chunked

class Prot:
    # TODO: Chain ID set to not none is not implemented
    # TODO: Does not handle multiple chains yet
    # TODO: Make coords the array for the 
    """
    Class that allows you query chemical properties of a protein
    """
    def __init__(self, pdb_path, top_path=None, chain_id=None):
        # Parse PDB file
        parser = PDBParser()
        structure = parser.get_structure("mol", pdb_path)
        # Grab first chain of first model as the protein
        if chain_id is None:
            self.chain = structure[0]["A"]
        # Grab initial Cartesian coordinates
        self.coords = []
        for atom in self.chain.get_atoms():
            atom_coords = atom.get_coord()
            self.coords.append(atom_coords)
        self.coords = np.array(self.coords)
        # Compute internal coordinates
        self.chain.atom_to_internal_coordinates()
        angle_types = ["phi", "psi", "chi1", "chi2", "chi3", "chi4", "chi5"]
        self.torsion_angles = []
        # Get torsion ids, initial torsion angle values, and if torsion is a chi angle
        for res in self.chain:
            if not res.internal_coord:
                print("WARNING: ", res, "could not be converted to internal coordinates")
            for angle_type in angle_types:
                angle = res.internal_coord.get_angle(angle_type)
                torsion_id = res.internal_coord.pick_angle(angle_type)
                is_chi = "chi" in angle_type
                self.torsion_angles.append([torsion_id, angle, is_chi])
        self.torsion_angles = np.array(self.torsion_angles)
        # Flags tell us if we need to update coords
        self.need_to_update_cartcoords = False
        self.need_to_update_intcoords = False
        return
    
    def get_charge(self, index: int):
        '''
        Returns the partial of the atom at a given index

        Args:
            index: the index of the atom (zero-indexed)
        '''
        return 

    
    def get_radius(self, index: int):
        '''
        Returns the radius of the atom at a given index

        Args:
            index: the index of the atom (zero-indexed)
        '''
        return
      
    def get_coords(self):
        '''
        Returns the Cartesian coordinates of the atom at a given index

        Args:
            index: the index of the atom (zero-indexed)
        '''
        if self.need_to_update_cartcoords:
            self.update_cartcoords_from_intcoords()
            self.need_to_update_cartcoords = False
        return self.coords

    def update_cartcoords_from_intcoords(self):
        """
        Updates the Cartesian coordinates to match the internal coordinates
        """
        self.chain.internal_to_atom_coordinates()
        for index, atom in enumerate(self.chain.get_atoms()):
            atom_coords = atom.get_coord()
            self.coords[index][0] = atom_coords[0]
            self.coords[index][1] = atom_coords[1]
            self.coords[index][2] = atom_coords[2]
        return
    
    def update_intcoords_from_cartcoords(self):
        """
        Updates internal coordinates to match the Cartesian coordinates
        """
        self.chain.atom_to_internal_coordinates()
        return

# TODO: Move to utils
def parse_biopython_torsion_id(torsion_id):
    """
    Parses biopythons torsion id.

    Returns:
        output: a 4 x 3 list containg the 4 atom ids (i.e. residue number, residue abbreviation, atom name)
    """
    atom_ids = torsion_id.split(":")
    output = []
    for index, atom_id in enumerate(atom_ids):
        atom_id_split = atom_id.split("_")
        output.append([atom_id_split[0], atom_id_split[1], atom_id_split[2]])
    return output
    
if __name__ == "__main__":
    pdb_path = "/home/conradli/SAC-for-H-Bond-Learning/data/pdb/1bdd.pdb"
    top_path = "/home/conradli/SAC-for-H-Bond-Learning/data/pdb/1bdd_proa.psf"
    protein = Prot(pdb_path, top_path)