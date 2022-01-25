import numpy as np

from collections.abc import Iterable
from mol_processors.Protein import Prot
from mol_processors.angle_lib import RotamerLibrary, BackboneDihedralLibrary

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
        torsion_indices = protein.select_torsions(res_index=index, selection="sidechain")
        protein.set_torsion_angles(torsion_indices, rot_sample)
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

    
    
    
