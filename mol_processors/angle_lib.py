import numpy as np
import os

from constants import RAMA_LINES_PER_NEIGHBOR, RAMA_LINES_PER_RESIDUE, RAMA_CHARS_PER_LINE, RAMA_RES_ORDER, RAMA_NEIGHBOR_ORDER

# TODO: non_rot_dist_type set to "density" and store_type set to "ram" not supported right now
class RotamerLibrary:
    """
    This class parses a rotamer library file (formatted like Dunbrack's) and
    provides functions to sample from the distributuion of rotamers
    """
    def __init__(self, rot_lib_dir, non_rot_dist_type="rotamer", store_type="disk"):
        """
        Creates a RotamerLibrary object given the directory with the rotamer library files

        Args:
            rot_lib_dir: the directory containing the rotamer library file
            non_rot_dist_type: the type of distribution to use for non-rotameric torsions (i.e. "rotamer" or "density")
            store_type: the type of storage to use for the rotamer library (i.e. "ram" or "disk")
        """
        self.non_rot_dist_type = non_rot_dist_type
        self.store_type = store_type
        self.rot_lib_dir = rot_lib_dir
        # Load the rotamer library into memory if storing in RAM
        if self.store_type == "ram":
            self.rot_lib = {}
            for filename in os.listdir(rot_lib_dir):
                # Only look at library files and rotamer libs (i.e. separated with three periods)
                filename_segs = filename.split('.')
                if filename_segs[-1] != 'lib' or len(filename_segs) != 4:
                    continue
                aa, _, rotlib_type, file = tuple(filename_segs)
 
                file = os.path.join(rot_lib_dir, filename)
                # Checking if it is a file
                if not os.path.isfile(file):
                    continue
        return
 
    def get_rotamer_dist(self, aa: str, phi: float, psi: float):
        """
        Returns possible rotamers and their probabilities given the phi and psi angles of a residue

        Args:
            aa: 3-letter amino acid abbreviation
            phi: phi angle in degrees
            psi: psi angle in degrees
        Returns:
            chi_probs: a ndarray of probabilities for sampling each rotamer for `aa`
            chi_samples: a ndarray of chi angles for each rotamer
        """
        # Dunbrack rotamers are in the format <aa>.bbdep.rotamers.lib
        aa = aa.lower()
        aa_rot_lib_path = os.path.join(self.rot_lib_dir, aa + ".bbdep.rotamers" + ".lib")
        if not os.path.isfile(aa_rot_lib_path):
            print("Could not find rotamer library for amino acid " + aa)
        with open(aa_rot_lib_path, "r") as aa_rot_lib:
            aa = aa.upper()
            # Round phi and psi to the nearest 10s digit
            psi = round(psi/10)*10
            phi = round(phi/10)*10
            chi_samples = []
            chi_probs = []
            chi_prob_sum = 0
            found_angles = False
            for line in aa_rot_lib:
                # Skip comments
                if line[0] == '#':
                    continue
                rot_entry = line.split()
                # Check if the phi and psi are in the rotamer library
                if rot_entry[1] == str(phi) and rot_entry[2] == str(psi):
                    # Get the frequency of the rotamer
                    chi_samples.append([float(rot_entry[i]) for i in range(9, len(rot_entry))])
                    chi_prob = float(rot_entry[8])
                    chi_probs.append(chi_prob)
                    chi_prob_sum += chi_prob
                    found_angles = True
                # Stop searching when we pass the psi/phi entry
                elif found_angles:
                    break
        # Normalize in case the sum is not exactly 1
        chi_probs = np.array(chi_probs)/chi_prob_sum
        return chi_probs, np.array(chi_samples)

# Disk not implemented
class BackboneDihedralLibrary:
    """
    This class parses a Dihedral library file (formatted like Dunbrack's) and
    provides functions to sample from the distributuion of backbone dihedrals
    """
    def __init__(self, backbone_lib_path: str):
        """
        Creates a BackboneDihedralLibrary object given the directory with the Ramachandran library files

        Args:
            backbone_lib_dir: the directory containing the Ramachandran library file
        """
        self.rama_left_path = os.path.join(backbone_lib_path, "rama_left.txt" )
        self.rama_right_path = os.path.join(backbone_lib_path, "rama_right.txt")
        if not os.path.isfile(self.rama_left_path) or not os.path.isfile(self.rama_right_path):
            raise FileNotFoundError("Could not find backbone dihedral library files")
        return
    
    def get_backbone_dist(self, aa: str, left_neighbor, right_neighbor):
        """
        Returns backbone dihedrals and their probabilities given 3 consecutive residues

        Args:
            aa: 3-letter abbreviation of central residue
            left_neighbor: 3-letter abbreviation of left residue
            right_neighbor: 3-letter abbreviation of right residue
        
        Returns:
            phi_psi_probs: a ndarray of probabilities of the ith phi-psi angles in `phi_psi_samples`
            phi_psi_samples: a N x 2 ndarray of phi-psi angle samples
        """
        rama_left = open(self.rama_left_path, "r")
        rama_right = open(self.rama_right_path, "r")
        # Seek in rama left and right to the first line for each residue
        aa_offset = RAMA_RES_ORDER[aa] * RAMA_LINES_PER_RESIDUE * RAMA_CHARS_PER_LINE
        left_neighbor_offset = aa_offset + (RAMA_NEIGHBOR_ORDER[left_neighbor] + 1) * RAMA_LINES_PER_NEIGHBOR * RAMA_CHARS_PER_LINE
        right_neighbor_offset = aa_offset + (RAMA_NEIGHBOR_ORDER[right_neighbor] + 1) * RAMA_LINES_PER_NEIGHBOR * RAMA_CHARS_PER_LINE
        # Grab left neighbor info
        rama_left.seek(left_neighbor_offset)
        left_info = np.zeros(RAMA_LINES_PER_NEIGHBOR)
        self.parse_rama_file(rama_left, left_neighbor, left_info) 
        # Grab right ALL info
        rama_right.seek(aa_offset)
        all_right_info = np.zeros(RAMA_LINES_PER_NEIGHBOR)
        self.parse_rama_file(rama_right, "ALL", all_right_info)
        # Grab right neighbor info
        rama_right.seek(right_neighbor_offset)
        right_info = np.zeros(RAMA_LINES_PER_NEIGHBOR)
        self.parse_rama_file(rama_right, right_neighbor, right_info)    
        # Calculate phi-psi probs conditioned on the three residues formula from Dunbrack
        p_star_sum = 0
        prob_distribution = np.zeros(RAMA_LINES_PER_NEIGHBOR)
        for index in range(RAMA_LINES_PER_NEIGHBOR):
            log_p_star_clr = -(left_info[index] + right_info[index] - all_right_info[index])
            p_star = np.exp(log_p_star_clr)
            prob_distribution[index] = p_star
            p_star_sum += p_star
        rama_left.close()
        rama_right.close()
        return prob_distribution / p_star_sum
        
    
    def parse_rama_file(self, rama_file, target_neighbor, output_array):
        """
        Helper function to parse a ramachandran file for a specific neighbor and store the results in `output_array`
        """
        # This for loop should never execute more than RAMA_LINES_PER_NEIGHBOR times
        for index, line in enumerate(rama_file):
            bb_angle_entry = line.split()
            neighbor = bb_angle_entry[2]
            if neighbor != target_neighbor:
                break
            # Store log_prob
            output_array[index] = float(bb_angle_entry[6])