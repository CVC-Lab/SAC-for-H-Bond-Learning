from typing import Iterable
from xmlrpc.client import Boolean
import MDAnalysis as mda
import numpy as np
from Bio.PDB import PDBParser
from constants import TORSION_ID_MAP, TORSION_TYPES
import numpy.typing as npt

class Prot:
    # TODO: Chain ID set to not none is not implemented
    # TODO: Does not handle multiple chains yet
    # TODO: Make coords the array for the 
    """
    Class that allows you query chemical properties of a protein 
    """
    coords: Iterable[float]
    """Atomic coordinates"""
    charges: Iterable[float]
    """Charges"""
    radii: Iterable[float]
    """Radii"""
    torsion_angles: Iterable[float]
    """Torsion angle values"""
    torsion_ids: Iterable[str]
    """Torsion ids"""
    torsion_mask: Iterable[int]
    """Torsion_ids has bunch of Nones where no angles exist (e.g. some amino acids do not have chi5). This masks those"""
    need_to_update_cartcoords: bool
    """Indicates if the internal coordinates changed and we need to change the cartesian coordinates"""
    need_to_update_intcoords: bool
    """Indicates if the Cartesian coordinates changed and we need to change the intenal coordinates"""

    def __init__(self, pdb_path, top_path=None, chain_id=None):
        # Parse PDB file
        parser = PDBParser()
        structure = parser.get_structure("mol", pdb_path)
        # Grab first chain of first model as the protein
        if chain_id is None:
            self.chain = structure[0]["A"]
        # Grab initial Cartesian coordinates, charges, and radii
        self.coords = []
        self.charges = []
        self.radii = []
        self.total_res = len(self.chain)
        for atom in self.chain.get_atoms():
            self.coords.append(atom.get_coord())
            self.charges.append(atom.get_charge())
            self.radii.append(atom.get_radius())
        self.coords = np.array(self.coords)
        self.charges = np.array(self.charges)
        self.radii = np.array(self.radii)
        # Compute internal coordinates
        self.chain.atom_to_internal_coordinates()
        self.torsion_ids = []
        self.torsion_angles = []
        # Store torsion ids, initial torsion angle values, and if torsion is a chi angle
        for res in self.chain:
            if not res.internal_coord:
                print("WARNING: ", res, "could not be converted to internal coordinates")
            for angle_type in TORSION_TYPES:
                angle = res.internal_coord.get_angle(angle_type)
                torsion_id = res.internal_coord.pick_angle(angle_type)
                self.torsion_ids.append(torsion_id)
                self.torsion_angles.append(angle)
        self.torsion_angles = np.array(self.torsion_angles)
        self.torsion_ids = np.array(self.torsion_ids)
        self.torsion_mask = self.torsion_ids != None
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
        return self.charges[index]

    
    def get_radius(self, index: int):
        '''
        Returns the radius of the atom at a given index

        Args:
            index: the index of the atom (zero-indexed)
        '''
        return self.radius[index]
      
    def get_coords(self):
        '''
        Returns the Cartesian coordinates of all the atoms
        '''
        if self.need_to_update_cartcoords:
            self.update_cartcoords_from_intcoords()
            self.need_to_update_cartcoords = False
        return self.coords

    def update_cartcoords_from_intcoords(self):
        """
        Updates the Cartesian coordinates to match the internal coordinates
        """
        # Update atomic coordinates
        self.chain.internal_to_atom_coordinates()
        # Update copy of atomic coordinates
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