from typing import Iterable
from xmlrpc.client import Boolean
import MDAnalysis as mda
import numpy as np
from Bio.PDB import PDBParser
from mol_processors.constants import TORSION_TYPES

class Prot:
    # TODO: Chain ID set to not none is not implemented
    # TODO: Does not handle multiple chains yet
    # TODO: Make coords the array for the 
    # TODO: Never modify the torsion angles array directly without calling set_torsion_angles
    """
    Class that allows you query chemical properties of a protein 
    """
    total_res: int
    """Total number of residues in chain"""
    residues: Iterable
    """Array of residues"""
    coords: Iterable[float]
    """Atomic coordinates"""
    charges: Iterable[float]
    """Charges"""
    radii: Iterable[float]
    """Radii"""
    torsion_angles: Iterable[float]
    """Torsion angle values. It will always be 7 * total_res (i.e. 7 possible torsion angles per residue)."""
    torsion_mask: Iterable[int]
    """torsion_angles has bunch of Nones where no angles exist (e.g. some amino acids do not have chi5). This masks those"""
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
        self.residues = list(self.chain.get_residues())
        for atom in self.chain.get_atoms():
            self.coords.append(atom.get_coord())
            self.charges.append(atom.get_charge())
            self.radii.append(atom.get_radius())
        self.coords = np.array(self.coords)
        self.charges = np.array(self.charges)
        self.radii = np.array(self.radii)
        # Compute internal coordinates
        self.chain.atom_to_internal_coordinates()
        self.torsion_angles = []
        # Store torsion ids, initial torsion angle values, and if torsion is a chi angle
        for res in self.chain:
            if not res.internal_coord:
                print("WARNING: ", res, "could not be converted to internal coordinates")
            for angle_type in TORSION_TYPES:
                angle = res.internal_coord.get_angle(angle_type)
                torsion_id = res.internal_coord.pick_angle(angle_type)
                self.torsion_angles.append(angle)
        self.torsion_angles = np.array(self.torsion_angles)
        self.torsion_mask = self.torsion_angles != None
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
        return self.coords
    
    def get_torsion_angles(self) -> np.ndarray:
        """
        Returns in-place torsion angles of the protein
        """
        if self.need_to_update_intcoords:
            self.update_intcoords_from_cartcoords()
        return self.torsion_angles
    
    def set_torsion_angles(self, indices:Iterable[int], new_angles: Iterable[float]) -> None:
        '''
        Sets the torsion angles given a set of residue numbers and angle keys but does not update
        the Cartesian coordinates. Must call `update_cartcoords_from_intcoords` to update cart coords.

        Args:
            protein: Protein chain we are altering
            indices: indices into the torsion ids that we want to change
            new_angles: the new angles for the torsion ids at `indices`
        Requirements:
            len(new_angles) = len(indices)
        '''
        # Update internal coordinates before modifying torsion angles
        if self.need_to_update_intcoords:
            self.update_intcoords_from_cartcoords()
            self.need_to_update_intcoords = False

        if len(indices) != len(new_angles):
            raise Exception("The length of new_angles and indices must be the same")
        # Set angles in biopython chain
        for index, torsion_index in enumerate(indices):
            active_res = self.residues[torsion_index // 7].internal_coord
            angle_key = TORSION_TYPES[torsion_index % 7]
            active_res.set_angle(angle_key, new_angles[index])  
        # Set angles in our copy of torsion angles
        self.need_to_update_cartcoords = True
        self.torsion_angles[indices] = np.copy(new_angles)

    def select_torsions(self, res_index: int = None, selection: str ="all"):
        """
        Returns torsion indices based whether we want just "backbone", "sidechains", or "all".
        
        Args:
            res_index: if set, restrict indices to a specific residue
            selection: 
                backbone: Only select phi and psi angles
                sidechains: Only select chi angles
                all: Select all angles

        Returns: an NDarray of torsion indices
        """
        total_torsions = self.total_res * 7
        indices = []
        # Return all indices
        if selection == "all":
            for i in range(total_torsions):
                if not self.torsion_mask[i]:
                        continue
                indices.append(i)   
            return indices
        # Get indices for specific residue
        if res_index != None:
            start_index = 7 * res_index
            if selection == "all":
                for i in range(start_index, start_index + 7):
                    if not self.torsion_mask[i]:
                        continue
                    indices.append(i)
                return indices
            elif selection == "backbone":
                for i in range(start_index, start_index + 2):
                    if not self.torsion_mask[i]:
                        continue
                    indices.append(i)
                return indices
            else:
                for i in range(start_index + 2, start_index + 7):
                    if not self.torsion_mask[i]:
                        continue
                    indices.append(i)
                return indices
        # Get indices out of whole chain
        else:
            for i in range(total_torsions):
                if not self.torsion_mask[i]:
                    continue
                elif selection == "backbone" and i % 7 < 2:
                    indices.append(i)
                elif selection == "sidechain" and i % 7 >= 2:
                    indices.append(i)
        return indices

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
        self.need_to_update_cartcoords = False
        return
    
    def update_intcoords_from_cartcoords(self):
        """
        Updates internal coordinates to match the Cartesian coordinates
        """
        self.chain.atom_to_internal_coordinates()
        self.need_to_update_intcoords = False
        return
    
if __name__ == "__main__":
    pdb_path = "/home/conradli/SAC-for-H-Bond-Learning/data/pdb/1bdd.pdb"
    top_path = "/home/conradli/SAC-for-H-Bond-Learning/data/pdb/1bdd_proa.psf"
    protein = Prot(pdb_path, top_path)