import numpy as np

from collections.abc import Iterable
from Protein import Prot
from angle_lib import RotamerLibrary, BackboneDihedralLibrary
from constants import TORSION_NUM_CHIS



def select_torsions(total_res, torsion_mask, res_index=None, selection="all"):
    """
    Returns torsion indices based whether we want just "backbone", "sidechains", or "all".
    
    Args:
        total_res: total number of residues in protein chain
        torsion_mask: indicates which elements of torsion id array do not contain torsion ids
        res_index: if set, restrict indices to a specific residue
        selection: 
            backbone: Only select phi and psi angles
            sidechains: Only select chi angles
            all: Select all angles

    Returns: an NDarray of torsion indices
    """
    total_torsions = total_res * 7
    # Return all indices
    if selection == "all":
        return np.arange(total_torsions)
    indices = []
    # Get indices for specific residue
    if res_index != None:
        start_index = 7 * int(selection)
        if selection == "all":
            for i in range(start_index, start_index + 7):
                if not torsion_mask[i]:
                    continue
                indices.append(i)
            return indices
        elif selection == "backbone":
            for i in range(start_index, start_index + 2):
                if not torsion_mask[i]:
                    continue
                indices.append(i)
            return indices
        else:
            for i in range(start_index + 2, start_index + 7):
                if not torsion_mask[i]:
                    continue
                indices.append(i)
            return indices
    # Get indices out of whole chain
    else:
        for i in range(total_torsions):
            if not torsion_mask[i]:
                continue
            elif selection == "backbone" and i % 7 < 2:
                indices.append(i)
            elif selection == "sidechain" and i % 7 >= 2:
                indices.append(i)
    return indices
            
def get_torsion_angle(protein: Prot, index: int) -> float:
    '''
    Returns the torsion angle for a given residue number and angle
    key

    Args:
        protein: the Protein object
        index: index of torsion angle
    
    Returns:
        a float representing the angle in radians 

    Runtime:
        O(1)
    '''
    return protein.torsion_angles[index] 

def set_torsion_angles(protein: Prot, indices:Iterable[int], new_angles: Iterable[float]) -> None:
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
    if len(indices) != len(new_angles):
        raise Exception("The length of new_angles and indices must be the same")
    torsion_ids = protein.torsion_ids[indices]
    for index, id in enumerate(torsion_ids):
        protein.chain.set_angle(id, new_angles[index])
    protein.need_to_update_cartcoords = True
    protein.torsion_angles[indices] = new_angles
    return

# TODO: Implement actually editing the pdb structure. Must be able to distinguish chi angles from phi and psi
def sample_dunbrack_sidechain_torsions(protein: Prot, rot_lib: RotamerLibrary, fixed_angles: Iterable[(int, str)] = []) -> None:
    for index, res in enumerate(protein.chain.get_residues()):
        res_name = res.get_resname().upper()
        phi = res.internal_coords.get_angle("phi")
        psi = res.internal_coords.get_angle("psi")
        # Sample from the posterior distribution
        rot_probs, rot_samples = rot_lib.get_rotamer_dist(res_name, phi, psi)
        rot_index = np.random.choice(len(rot_samples), p=rot_probs)
        rot_gaussian = rot_samples[rot_index]
        # Sample chi values from the Gaussian 
        rot_mean = rot_gaussian[:4]
        rot_std = rot_gaussian[4:]
        rot_sample = np.random.normal(rot_mean, rot_std)
        torsion_indices = select_torsions(protein.total_res, protein.torsion_mask, res_index=index, selection="sidechain")
        set_torsion_angles(protein, torsion_indices, rot_sample)
    return

# TODO: Implement actually editing the pdb structure
def sample_dunbrack_backbone_dihedrals(protein: Prot, dihedral_lib: BackboneDihedralLibrary, fixed_res: Iterable[int] = []) -> None:
    phi_psi_probs = dihedral_lib.get_backbone_dist("arg", "val", "leu")
    # Sample the lower left corner of 5 x 5 phi-psi grid
    phi_psi_index = np.random.choice(len(phi_psi_probs), p=phi_psi_probs)
    lower_left_phi = -180 + 5 * (phi_psi_index // 72) # phi increases by 5 degrees every 72 lines
    lower_left_psi = -180 + 5 * (phi_psi_index % 72) # psi increases by 5 degrees but wraps around every 72 lines
    # Sample uniformly from 5 x 5 grid with phi_psi_sample at lower left corner
    phi_sample = np.random.uniform(lower_left_phi, lower_left_phi + 5)
    psi_sample = np.random.uniform(lower_left_psi, lower_left_psi + 5)
    return

if __name__ == "__main__":
    import tracemalloc
    tracemalloc.start()
    rot_lib = "/home/conradli/SAC-for-H-Bond-Learning/data/torsion_libs/Dunbrack2011-5/ExtendedOpt1-5"
    bb_lib = "/home/conradli/SAC-for-H-Bond-Learning/data/torsion_libs/ramachandran/"
    aa = "HIS"
    lib = BackboneDihedralLibrary(bb_lib)

    
    
    
