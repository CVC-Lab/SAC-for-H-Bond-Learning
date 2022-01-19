import MDAnalysis as mda
from Bio.PDB import PDBParser

class Prot:
    # TODO: Chain ID set to not none is not implemented
    # TODO: Does not handle multiple chains yet
    """
    Class that allows you query chemical properties of a protein
    """
    def __init__(self, pdb_path, top_path=None, chain_id=None):
        # Parse PDB file
        parser = PDBParser()
        structure = parser.get_structure("mol", pdb_path)
        # Grab first chain of first model as the protein
        if chain_id is None:
            chain = structure[0][0]
        # Grab atoms from chain
        self.atoms = chain.get_atoms()
        # Flags tell us if we need to update coords
        self.need_to_update_cartcoords = False
        self.need_to_update_intcoords = False
        self.fixed_angles = []
        return
    
    def get_charge(self, index: int):
        '''
        Returns the partial of the atom at a given index

        Args:
            index: the index of the atom (zero-indexed)
            
        '''
        return self.atoms[index].get_charge()

    
    def get_radius(self, index: int):
        '''
        Returns the radius of the atom at a given index

        Args:
            index: the index of the atom (zero-indexed)
        '''
        return self.atoms[index].get_radius()
    
    
    def get_coords(self, index: int):
        '''
        Returns the Cartesian coordinates of the atom at a given index

        Args:
            index: the index of the atom (zero-indexed)
        '''
        return self.atoms[index].get_coords()

    # Updates 
    def update_cartcoords_from_intcoords(self):
        return
    
    def update_intcoords_from_cartcoords(self):
        return
    
    def get_coords(self, pdb_path, top_path=None):
        if self.need_to_update_cartcoords:
            self.update_cartcoords_from_intcoords()
        return self.coords
    