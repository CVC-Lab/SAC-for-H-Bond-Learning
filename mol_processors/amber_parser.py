import numpy as np

class PrmtopParser:
    """
    charges (NUM_ATOMS): charge of each atom
    masses (NUM_ATOMS): mass of each atom
    atom_names (NUM_ATOMS): name of each atom
    atom_numbers (NUM_ATOMS): atomic number of each atom
    atomic_radii (NUM_ATOMS): atomic radii
    bonds (NUM_BONDS x 3): Bonds without hydrogen come first. Each bonds is defined by Atom 1, Atom 2, Bond Parameter Index
    num_bonds_with_hydro (int): Number of bonds with hydrogen
    angles (NUM_ANGLES x 3): Bonds without hydrogen come first. Each angle is defined by Atom 1, Atom 2, Atom 3, Angle Parameter Index
    num_angles_with_hydro (int): Number of angles with hydrogen
    dihedrals (NUM_DIHEDRALS x 4): Atom 1, Atom 2, Atom 3, Atom 4, Dihedral Parameter Index
    num_dihedrals_with_hydro (int): Number of dihedrals with hydrogen
    lj_params (NUM_UNIQUE_ATOM_TYPE_PAIRS x 2): Lennard-Jones A and B parameters
    lj_pair_index (NUM_UNIQUE_ATOM_TYPE_PAIRS): NONBONDED_PARM_INDEX[NTYPES x (ATOM_TYPE_INDEX(i) - 1) + ATOM_TYPE_INDEX(j)
    lj_atom_type_index (NUM_LJ_ATOM_TYPES): LJ Atom type index used with lj pair index to find LJ params
    excluded_atoms (NUM_ATOMS * NUM_EXCLUDED_ATOMS[ATOM]): list of excluded atoms
    num_excluded_atoms (NUM_ATOMS): Number of atoms that are excluded from non-bonded interactions for each atom
    """
    def __init__(self, prmtop_path):
        prmtop_dict = self.parse_prmtop(prmtop_path)
        # Atomic properties
        self.charges, self.masses, self.atom_names, self.atom_numbers, self.atomic_radii \
            = self.get_atom_properties(prmtop_dict)
        # Bonds, Angles, and Dihedrals
        self.bonds, self.num_bonds_with_hydro = self.get_bonds(prmtop_dict)
        self.angles, self.num_angles_with_hydro = self.get_angles(prmtop_dict)
        self.dihedrals, self.num_dihedrals_with_hydro = self.get_dihedrals(prmtop_dict)
        # Bond force field params
        self.bond_force, self.bond_equil = self.get_bond_length_params(prmtop_dict)
        # Angle force field params
        self.angle_force, self.angle_equil = self.get_bond_angle_params(prmtop_dict)
        # dihedral force field params
        self.dihedral_force, self.dihedral_periodicity, self.dihedral_phase \
            = self.get_dihedral_angle_params(prmtop_dict)
        # LJ force field params
        self.lj_params = self.get_lj_params(prmtop_dict)
        self.lj_atom_type_index, self.lj_pair_index = self.get_lj_indices(prmtop_dict)
        # Excluded atoms and number of excluded atoms for each atom
        self.excluded_atoms, self.num_excluded_atoms= self.get_excluded_atoms(prmtop_dict)
        del prmtop_dict
        return

    def get_atom_properties(self, prmtop_dict):
        """
        Retrieves charges, masses, atom names, atom numbers, and atomic radii
        """
        charges = np.array(prmtop_dict["charge"], dtype=np.float32)
        masses = np.array(prmtop_dict["mass"], dtype=np.float32)
        atom_names = prmtop_dict["atom_name"]
        atom_numbers = np.array(prmtop_dict["atomic_number"], dtype=np.int32)
        atomic_radii = np.array(prmtop_dict["radii"], dtype=np.float32)
        return charges, masses, atom_names, atom_numbers, atomic_radii

    def get_bonds(self, prmtop_dict):
        """
        Retrieves bonds with hydrogen and bonds without hydrogen from prmtop dictionary
        Extracts the two atoms and the bond parameter index associated with each bond
        """
        # Validate bond arrays
        if len(prmtop_dict["bonds_inc_hydrogen"]) % 3 != 0:
            raise Exception("Bonds_inc_hydrogen must have a length divisible by 3")
        if len(prmtop_dict["bonds_without_hydrogen"]) % 3 != 0:
            raise Exception("Bonds_without_hydrogen must have a length divisible by 3")
        # Load bonds into numpy arrays
        bonds_with_hydro = np.array(prmtop_dict["bonds_inc_hydrogen"], dtype=np.int32)
        bonds_without_hydro = np.array(prmtop_dict["bonds_without_hydrogen"], dtype=np.int32)
        num_bonds_with_hydro = len(bonds_with_hydro) // 3
        bonds = np.concatenate((bonds_with_hydro, bonds_without_hydro))
        # Reshape bond array so each row is a bond
        bonds = np.reshape(bonds, (len(bonds) // 3, 3))
        # Correct indices to match true atomic indices (A = N/3 + 1)
        for bond in bonds:
            assert bond[0] % 3 == 0
            bond[0] = bond[0] // 3 + 1
            assert bond[1] % 3 == 0
            bond[1] = bond[1] // 3 + 1
        return bonds, num_bonds_with_hydro
    
    def get_angles(self, prmtop_dict):
        """
        Extracts bond angles from prmtop dictionary and corrects atomic coordinates
        """
        # Validate bond arrays
        if len(prmtop_dict["angles_inc_hydrogen"]) % 4 != 0:
            raise Exception("Angles_inc_hydrogen must have a length divisible by 4")
        if len(prmtop_dict["angles_without_hydrogen"]) % 4 != 0:
            raise Exception("Angles_without_hydrogen must have a length divisible by 4")
        # Load angles into numpy arrays
        angles_with_hydro = np.array(prmtop_dict["angles_inc_hydrogen"], dtype=np.int32)
        angles_without_hydro = np.array(prmtop_dict["angles_without_hydrogen"], dtype=np.int32)
        num_angles_with_hydro = len(angles_with_hydro) // 4
        angles = np.concatenate((angles_with_hydro, angles_without_hydro))
        # Reshape angle array so each row is a bond
        angles = np.reshape(angles, (len(angles) // 4, 4))
        # Correct indices to match true atomic indices (A = N/3 + 1)
        for angle in angles:
            assert angle[0] % 3 == 0
            angle[0] = angle[0] // 3 + 1
            assert angle[1] % 3 == 0
            angle[1] = angle[1] // 3 + 1
            assert angle[2] % 3 == 0
            angle[2] = angle[2] // 3 + 1
        return angles, num_angles_with_hydro

    def get_dihedrals(self, prmtop_dict):
        """
        Extracts dihedral angles from prmtop dictionary and corrects atomic coordinates
        """
        # Validate bond arrays
        if len(prmtop_dict["dihedrals_inc_hydrogen"]) % 5 != 0:
            raise Exception("dihedrals_inc_hydrogen must have a length divisible by 5")
        if len(prmtop_dict["dihedrals_without_hydrogen"]) % 5 != 0:
            raise Exception("dihedrals_without_hydrogen must have a length divisible by 5")
        # Load dihedrals into numpy arrays
        dihedrals_with_hydro = np.array(prmtop_dict["dihedrals_inc_hydrogen"], dtype=np.int32)
        dihedrals_without_hydro = np.array(prmtop_dict["dihedrals_without_hydrogen"], dtype=np.int32)
        num_dihedrals_with_hydro = len(dihedrals_with_hydro) // 5
        dihedrals = np.concatenate((dihedrals_with_hydro, dihedrals_without_hydro))
        # Reshape angle array so each row is a bond
        dihedrals = np.reshape(dihedrals, (len(dihedrals) // 5, 5))
        # Correct atom indices to match true atomic indices (A = N/3 + 1)
        for dihedral in dihedrals:
            assert dihedral[0] % 3 == 0
            dihedral[0] = dihedral[0] // 3 + 1
            assert dihedral[1] % 3 == 0
            dihedral[1] = dihedral[1] // 3 + 1
            assert dihedral[2] % 3 == 0
            dihedral[2] = dihedral[2] // 3 + 1
            assert dihedral[3] % 3 == 0
            dihedral[3] = dihedral[3] // 3 + 1
        return dihedrals, num_dihedrals_with_hydro

    def get_lj_params(self, prmtop_dict):
        """
        Retrieves LJ parameters
        """
        a_coeff = np.array(prmtop_dict["lennard_jones_acoef"], dtype=np.float32)
        b_coeff = np.array(prmtop_dict["lennard_jones_bcoef"], dtype=np.float32)
        return np.vstack((a_coeff, b_coeff))
    
    def get_lj_indices(self, prmtop_dict):
        """
        Retrieves LJ atom type indices and LJ pair indices used find LJ parameters a
        given nonbonded pair of atoms
        """
        lj_pair_indices = np.array(prmtop_dict["nonbonded_parm_index"], dtype=np.int32)
        lj_atom_type_indices = np.array(prmtop_dict["atom_type_index"], dtype=np.int32)
        return lj_atom_type_indices, lj_pair_indices

    def get_excluded_atoms(self, prmtop_dict):
        """
        Retrieves indices of atoms that are excluded from certain non-bonded \
        interactions (e.g. 1-4 and 1-2)
        """
        excluded_atoms = np.array(prmtop_dict["excluded_atoms_list"], dtype=np.int32)
        num_excluded_atoms = np.array(prmtop_dict["number_excluded_atoms"], dtype=np.int32)
        return excluded_atoms, num_excluded_atoms

    def get_bond_length_params(self, prmtop_dict):
        """
        Retrieves bond length parameters
        """
        bond_force = np.array(prmtop_dict["bond_force_constant"], dtype=np.float32)
        bond_equil = np.array(prmtop_dict["bond_equil_value"], dtype=np.float32)
        return bond_force, bond_equil
    
    def get_bond_angle_params(self, prmtop_dict):
        """
        Retrieves bond angle parameters
        """
        angle_force = np.array(prmtop_dict["angle_force_constant"], dtype=np.float32)
        angle_equil = np.array(prmtop_dict["angle_equil_value"], dtype=np.float32)
        return angle_force, angle_equil
    
    def get_dihedral_angle_params(self, prmtop_dict):
        """
        Retrieves dihedral angle parameters
        """
        dihedral_force = np.array(prmtop_dict["dihedral_force_constant"], dtype=np.float32)
        dihedral_periodicity = np.array(prmtop_dict["dihedral_periodicity"], dtype=np.float32)
        dihedral_phase = np.array(prmtop_dict["dihedral_phase"], dtype=np.float32)
        return dihedral_force, dihedral_periodicity, dihedral_phase

    def parse_prmtop(self, prmtop_path):
        """
        Parse a prmtop file and returns a dictionary containing all flag information
        """
        prmtop_file = open(prmtop_path, "r")
        prmtop_lines = prmtop_file.readlines()
        prmtop_file.close()
        prmtop_dict = self.setup_prmtop_dict()
        current_key = None
        for line in prmtop_lines:
            # Check for new section
            if line.startswith("%"):
                if line.startswith("%FLAG"):
                    flag = line.split()[1].lower()
                    # Only add the sections we are interested in
                    if flag in prmtop_dict.keys():
                        current_key = flag
                        print("Key Set to", flag)
                    else:
                        print("Skipping flag: {} in prmtop file".format(flag))
                        current_key = None
                # Skip the format line
                elif line.startswith("%FORMAT"):
                    continue
            # Add info to list at the current key
            elif current_key != None:
                prmtop_dict[current_key].extend(line.split())
        return prmtop_dict
    
    def setup_prmtop_dict(self):
        prmtop_dict = {}
        # Set up prmtop dictionary with empty lists
        prmtop_dict["pointers"] = []
        prmtop_dict["atom_name"] = []
        prmtop_dict["charge"] = []
        prmtop_dict["atomic_number"] = []
        prmtop_dict["mass"] = []
        prmtop_dict["atom_type_index"] = []
        prmtop_dict["nonbonded_parm_index"] = []
        prmtop_dict["residue_label"] = []
        prmtop_dict["residue_pointer"] = []
        prmtop_dict["bond_force_constant"] = []
        prmtop_dict["bond_equil_value"] = []
        prmtop_dict["angle_force_constant"] = []
        prmtop_dict["angle_equil_value"] = []
        prmtop_dict["dihedral_force_constant"] = []
        prmtop_dict["dihedral_periodicity"] = []
        prmtop_dict["dihedral_phase"] = []
        prmtop_dict["scee_scale_factor"] = []
        prmtop_dict["scnb_scale_factor"] = []
        prmtop_dict["lennard_jones_acoef"] = []
        prmtop_dict["lennard_jones_bcoef"] = []
        prmtop_dict["bonds_inc_hydrogen"] = []
        prmtop_dict["bonds_without_hydrogen"] = []
        prmtop_dict["angles_inc_hydrogen"] = []
        prmtop_dict["angles_without_hydrogen"] = []
        prmtop_dict["dihedrals_inc_hydrogen"] = []
        prmtop_dict["dihedrals_without_hydrogen"] = []
        prmtop_dict["excluded_atoms_list"] = []
        prmtop_dict["number_excluded_atoms"] = []
        prmtop_dict["radii"] = []
        return prmtop_dict

def validate_prmtop_parser(prmtop_parser):
    """
    Checks to make sure that the prmtop parsed correctly
    """
    return

if __name__ == "__main__":
    prmtop_path = "/home/conradli/SAC-for-H-Bond-Learning/data/pdb/dialanine.top"
    prmtop_parser = PrmtopParser(prmtop_path)
