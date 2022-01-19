import numpy as np
import os

from collections.abc import Iterable
from Bio.PDB import internal_coords
from Protein import Prot
from pdb_processor import res_3to1

class RotamerLibrary:
    """
    This class parses a rotamer library file (formatted like Dunbrack's) and
    provides functions to sample from the distributuion of rotamers

    Args:
        rot_lib_dir: the directory containing the rotamer library file
        non_rot_dist_type: the type of distribution to use for non-rotameric torsions
    """
    def __init__(self, rot_lib_dir, non_rot_dist_type="rotamer"):
        for filename in os.listdir(rot_lib_dir):
            aa, _, lib_type, file_ending = filename.split('.')
            print(aa, lib_type, file_ending)
            # Skip if not rotamer library file
            if file_ending != 'lib':
                continue
            file = os.path.join(rot_lib_dir, filename)
            # Checking if it is a file
            if not os.path.isfile(file):
                continue
            
        return

class BackboneDihedralLibrary:
    """
    This class parses a Dihedral library file (formatted like Dunbrack's) and
    provides functions to sample from the distributuion of backbone dihedrals
    """
    def __init__(self):
        return

def get_torsion_angle(protein: Prot, res_num: int, angle_key: str) -> float:
    '''
    Returns the torsion angle for a given residue number and angle
    key

    Args:
        protein: the Protein object
        res_num: a 1-indexed residue number 
        angle_key: an angle key ("phi", "psi", "omega", "chi1", "chi2", "chi3", "chi4", "chi5)
    
    Returns:
        a float representing the angle in radians 

    Runtime:
        O(1)
    '''
    return 

def set_torsion_angles(protein: Prot, res_num: Iterable[int], angle_key: Iterable[str]) -> None:
    '''
    Sets the torsion angles given a set of residue numbers and angle keys but does not update
    the Cartesian coordinates. Must call `update_cartcoords_from_intcoords` to update cart coords.
    '''
    return

def sample_dunbrack_sidechain_torsions(protein: Prot, rot_lib: RotamerLibrary, fixed_angles: Iterable[(int, str)] = []) -> None:
    return

def sample_dunbrack_backbone_dihedrals(protein: Prot, dihedral_lib:  BackboneDihedralLibrary, fixed_res: Iterable[int] = []) -> None:
    return


if __name__ == "__main__":
    lib = RotamerLibrary("/home/conradli/SAC-for-H-Bond-Learning/data/torsion_libs/Dunbrack2011-5/ExtendedOpt1-5")