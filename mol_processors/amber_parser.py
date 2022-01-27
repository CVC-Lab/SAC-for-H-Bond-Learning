import numpy as np

class PrmtopParser:
    
    def __init__(self, prmtop_path):
        prmtop_dict = self.parse_prmtop(prmtop_path)
        # Atomic properties
        self.charges, self.masses, self.atom_names, self.atom_numbers, self.atomic_radii \
            = self.get_atom_properties(prmtop_dict)
        # Bonds, Angles, and Dihedrals
        self.bonds = None
        self.angles = None
        self.dihedrals = None
        # Force field parameters
        self.bond_params = None
        self.angle_params = None
        self.torsion_params = None
        self.lj_params = None
        self.atom_lj_index = None
        # Excluded atoms
        self.excluded_atoms = None
        self.num_excluded_atoms = None
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
    
    def get_nonbonded_params(self, prmtop_dict):
        """
        Retrieves nonbonded parameters
        """
        
        return a_coeff, b_coeff
    
    def get_bond_length_params(self, bond_index):
        """
        Retrieves bond length parameters
        """
        return bond_force, bond_equil
    
    def get_bond_angle_params(self, bond_angle_index):
        """
        Retrieves bond angle parameters
        """
        return angle_force, angle_equil
    
    def get_torsion_angle_params(self, torsion_index):
        """
        Retrieves torsion angle parameters
        """
        return torsion_force, torsion_periodicity, torsion_phase
    
    def parse_prmtop(self, prmtop_path):
        """
        Parse a prmtop file and returns a dictionary containing all flag information
        """
        prmtop_file = open(prmtop_path, "r")
        prmtop_lines = prmtop_file.readlines()
        prmtop_file.close()
        prmtop_dict = self.setup_prmtop_dict()
        for line in prmtop_lines:
            current_key = None
            # Check for new section
            if line.startswith("%"):
                if line.startswith("%FLAG"):
                    flag = line.split()[1]
                    print("Considering flag:", flag)
                    # Only add the sections we are interested in
                    if current_key in prmtop_dict.keys():
                        current_key = flag
                    else:
                        print("Skipping flag: {} in prmtop file".format(flag))
                        current_key = None
                # Skip the format line
                elif line.startswith("%FORMAT"):
                    continue
            # Add info to list at the current key
            elif current_key != None:
                prmtop_dict[current_key].extend(line.split())
    
    def setup_prmtop_dict(self):
        prmtop_dict = {}
        # Set up prmtop dictionary with empty lists
        prmtop_dict["pointers"] = []
        prmtop_dict["atom_name"] = []
        prmtop_dict["charge"] = []
        prmtop_dict["atomic_number"] = []
        prmtop_dict["mass"] = []
        prmtop_dict["atom_type_index"] = []
        prmtop_dict["number_excluded_atoms"] = []
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
        prmtop_dict["lennard_jones_acoeff"] = []
        prmtop_dict["lennard_jones_bcoeff"] = []
        prmtop_dict["bonds_inc_hydrogen"] = []
        prmtop_dict["bonds_without_hydrogen"] = []
        prmtop_dict["angles_inc_hydrogen"] = []
        prmtop_dict["angles_without_hydrogen"] = []
        prmtop_dict["dihedrals_inc_hydrogen"] = []
        prmtop_dict["dihedrals_without_hydrogen"] = []
        prmtop_dict["excluded_atoms_list"] = []
        prmtop_dict["radii"] = []
        return prmtop_dict
    

if __name__ == "__main__":
    prmtop_path = "/home/conradli/SAC-for-H-Bond-Learning/data/pdb/dialanine.top"
    prmtop_parser = PrmtopParser(prmtop_path)
