# Imports
import random
from utils import write_array
from Octree.octree import Protein
import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader
# Rosetta imports
from pyrosetta import *
from pyrosetta.toolbox import *
from rosetta.core.pack.task import TaskFactory
from rosetta.core.pack.task import operation
from pyrosetta.rosetta.core.scoring import ScoreType, Ramachandran, ScoringManager
from pyrosetta.rosetta.protocols.backbone_moves import RandomizeBBByRamaPrePro
from pyrosetta.rosetta.protocols.relax import FastRelax
import pyrosetta.rosetta.core.pack.rotamer_set as rotamer_set
from pyrosetta.rosetta.core.pack import create_packer_graph
from pyrosetta.rosetta.protocols import minimization_packing
from pyrosetta.rosetta.core.id import TorsionID, TorsionType
import pyrosetta.rosetta.protocols as protocols
#from .bat import BAT
import nglview as nv
import numpy as np

# Dictionary used for one hot encoding of element names
one_hot_map = {"C": 0, "N": 1, "H": 2, "O": 3, "S": 2}
NUM_ELEMENTS = 5

# Class that allows you to query data about a specific protein and edit its
# configuration
class Prot:

    # Initializes the Prot Object
    # pdb_file: path to pdb
    # psf_file: path to psf
    def __init__(self, pdb_file=None, psf_file=None, mol2_file=None, prm_file=None, rtf_file=None, aprm_file=None, pdb_id=None, seq=None, rot_lib_type="independent"):
        '''
        MDANALYSIS WAY
        # Load in pdb and psf
        self.prot = mda.Universe(psf_file, psf_file)
        protein_residues = self.prot.select_atoms("protein")
        # Store the Cartesian coordinates
        self.cart = None
        # Create a bond adjacency list
        # Calculates a N x 4 numpy matrix representing the internal coordinates
        #self.intern = BAT(protein_residues) # N x 4 matrix representing the internal coordinates 
        #self.intern.run()
        # Backbone atoms
        '''
        # Load in pdb file
        if pdb_file != None:
            print("Loading pdb fie into Rosetta...")
            self.pose = pose_from_pdb(pdb_file)
            self.pdb_name = pdb_file.split("/")[-1].split(".pdb")[0]
            # if we have the files create an octree
            if psf_file != None and mol2_file != None and prm_file != None and rtf_file != None and aprm_file != None:
                outnFile = self.pdb_name + "-outn.pdb"
                outoFile = self.pdb_name + "-outo.pdb"
                self.octree_rep = Protein(pdb_file, psf_file, mol2_file, prm_file, rtf_file, aprm_file, outnFile, outoFile)
        # Grab the pdb id straight from RCSB
        elif pdb_id != None:
            print("Loading pdb id into Rosetta")
            self.pose = pose_from_rcsb(pdb_id)
            self.pdb_name = pdb_id
        # Default to sequence
        elif seq != None:
            self.pose = pose_from_sequence(seq, auto_termini=True)
            self.pdb_name = "alanine_dipeptide"

        # Atom ID to pdb_atom_id dict
        self.incl_adj = False

        self.atom_id_to_pdb_id = {}
        # Store pdb file name and psf file name
        self.pdb_file = pdb_file
        self.psf_file = psf_file
        # Get total atoms and total residues
        self.num_residues = self.pose.total_residue()
        self.num_atoms = self.pose.total_atoms()
        self.num_virt_atoms = self.get_num_virt_atoms()
        self.ind_rot_set = None
        # Store backbone independent sets if needed
        if rot_lib_type == "independent":
            self.setup_ind_rot_set()
        # Setup rotamer packer task
        self.rot_set_fact = rotamer_set.RotamerSetFactory()
        self.pack_task = standard_packer_task(self.pose)
        self.pack_task.restrict_to_repacking() # Do not allow sequence changes
        # Store atom ids
        self.atom_ids = self.get_atom_ids()
        # Store secondary structure per atom
        self.update_sec_struct()
        # Store Cartesian coordinates (do not store virtual atoms) 
        self.cart_coords = np.zeros((self.num_atoms-self.num_virt_atoms, 3))
        # Store bonded adjacency list
        self.bond_adj_mat = self.get_bond_adj_mat()
        # Store hbond set. This only gets filled when self.get_hbonds is called
        self.hbond_set = None
        # Generates chemical features
        self.atom_chem_features = self.generate_chemical_features()
        self.update_cart_coords()
        # Set the scoring functon
        self.score_function = self.set_score_function(only_hbond=False)
        print("reached")
        self.score_function.setup_for_packing(self.pose, self.pack_task.repacking_residues(), self.pack_task.designing_residues())
        # Set up packer graph now that score function is set
        self.pack_graph = create_packer_graph(self.pose, self.score_function, self.pack_task)
        # Save original pose
        self.original_pose = self.pose.clone()
        

    ################# INPUT/OUTPUT FUNCTIONS #################
    
    # Writes the protein configuration to a PDB file
    def write_to_pdb(self, output_file):
        self.pose.dump_pdb(output_file)
        return
    
    # Writes numpy coords to pdb files
    # Coord_files is a list of coord files of saved npy files
    # NOT SURE IF THIS ACTUALLY WORKS YET
    def np_coords_to_pdbs(self, coord_files, output_dir, coord_file_type= "npy"):
        # Check for valid file type
        if coord_file_type != "npy" or coord_file_type != "text":
            raise Exception("Invalid file type")
        # Loop through all coord files
        for index, coord_file in enumerate(coord_files):
            coords = []
            if coord_file_type == "npy":
                coords = np.load(coord_file)
            elif coord_file_type == "text":
                coords = np.loadtxt(coord_file)
            # Convert each coord to xyz vector
            xyz = []
            for coord in coords:
                xyz.append(numeric.xyzVector_float(coord))
            # Set the coords
            self.pose.batch_set_xyz(Vector1(self.atom_ids), Vector1(xyz))
            output_file_name = self.pdb_name + str(index) + ".pdb"
            self.write_to_pdb(os.path.join(output_dir, output_file_name))
        return 

    # Coords can be a pdb file, traj file (dcd, xtc,...), or a numpy array (frame, num_atoms, 3)
    def visualize_protein(self, coords=None, default=None, default_representation=False):
        # Default is to use the cuurent Pose coords
        if coords == None:
            w = nv.show_rosetta(self.pose)
            return w
        # Load in np array of coords
        elif type(coords).__name__ == "ndarray":
            u = mda.Universe(self.psf_file, coords, format=MemoryReader, order='fac')
        # Load in traj or pdb file 
        elif type(coords).__name__ == "str":
            # Select all atoms associated with a protein
            u = mda.Universe(self.psf_file, coords)
        else:
            raise Exception("coords must be either a string or numpy array")
        protein_residues = u.select_atoms("protein")
        w = nv.show_mdanalysis(protein_residues, default=default, default_representation=default_representation)
        return w
    
    # Writes the Cartesian coordinates to disk
    def write_cart_coords(self, output_file, file_type=None):
        write_array(self.cart_coords, output_file, file_type)
        return        

    ################# ATOM IDS #################
    # Get the number of virtual atoms
    def get_num_virt_atoms(self):
        count = 0
        for res_index in range(1, self.num_residues + 1):
            residue = self.pose.residue(res_index)
            num_res_atoms = residue.natoms()
            # Loop through all atoms of each residue
            for res_atom_index in range(1, num_res_atoms+1):
                if(residue.atom_type(res_atom_index).is_virtual()):
                    count += 1
        return count

    # Returns a numpy array of atom ids in 
    # Atom ID stores the residue id and atom number: 
    # 1. atom_id.rsd()
    # 2. atom_id.atomno()
    def get_atom_ids(self):
        atom_ids = []
        pdb_atom_id = 0
        for res_index in range(1, self.num_residues + 1):
            residue = self.pose.residue(res_index)
            num_res_atoms = residue.natoms()
            # Loop through all atoms of each residue
            for res_atom_index in range(1, num_res_atoms+1):
                # skip virtual atoms (i.e. NV on proline)
                if(not residue.atom_type(res_atom_index).is_virtual()):
                    atom_ids.append(AtomID(res_atom_index, res_index))
                    self.atom_id_to_pdb_id[(res_atom_index, res_index)] = pdb_atom_id
                    pdb_atom_id += 1
        return atom_ids

    
    # Returns the internal coordinates
    # Converts the atomsitic positions to internal coordinates using mdanalysis
    # In interanl coordinates, each atom is defined by
    # 1. Bond Length (defined by 2 atoms)
    # 2. Bond Angle (defined by 3 atoms)
    # 3. Dihedral Angle (defined by 4 atoms)
    # TODO
    def get_internal_coords(self):
        return None
    
    ################# CARTESIAN COORDINATES #################

    # Returns the Cartesian coordinates
    def get_cart_coords(self):
        return self.cart_coords

    # Returns the Cartesian coordinates or the X, Y, Z position of each atom
    def update_cart_coords(self):
        pdb_atom_index = 0
        # Loop through all residues
        for atom_id in self.atom_ids:
                res_id = atom_id.rsd()
                atom_index = atom_id.atomno()
                xyz = self.pose.residue(res_id).atom(atom_index).xyz()
                # Update cart coords in cart coords
                self.cart_coords[pdb_atom_index, 0] = xyz[0]
                self.cart_coords[pdb_atom_index, 1] = xyz[1]
                self.cart_coords[pdb_atom_index, 2] = xyz[2]
                # Update cart coords in atom features
                self.atom_chem_features[pdb_atom_index, 0] = xyz[0]
                self.atom_chem_features[pdb_atom_index, 1] = xyz[1]
                self.atom_chem_features[pdb_atom_index, 2] = xyz[2]
                pdb_atom_index += 1
    
    ################# SECONDARY STRUCTURES #################

    # Gets the DSSP secondary structure assignment of the residue
    def get_sec_struct(self, residue_index):
        return self.pose.secstruct(residue_index)

    # Updates the scondary structure data
    def update_sec_struct(self):
        self.pose.display_secstruct()

    ################# COVALENT BONDS #################

    # The adjacency list showing which atoms are bonded to which atoms
    def get_bond_adj_list(self):
        bond_adj_list = {}
        # Loop through all residues
        for atom_id in self.atom_ids:
            neighbor_atom_ids  = self.pose.conformation().bonded_neighbor_all_res(atom_id)
            bond_adj_list[atom_id] = np.array(neighbor_atom_ids)
        return bond_adj_list
    
    # The adjacency matrix showing which atoms are bonded to which atoms
    def get_bond_adj_mat(self, self_neighbors=True):
        bond_adj_list = []
        # Loop through all residues
        for atom_id in self.atom_ids:
            bond_adj_list.append(self.get_bond_adj_for_atom(atom_id))
        if self_neighbors:
            return np.array(bond_adj_list) + np.eye(len(self.atom_ids), len(self.atom_ids))
        return np.array(bond_adj_list)

    # Can only be called after get_bond_adj_list
    def get_bond_adj_for_atom(self, atom_id):
        adj_info =  [0] * len(self.atom_ids)
        neighbors = self.pose.conformation().bonded_neighbor_all_res(atom_id)
        # Query for pdb atom id of neighbor
        for neighbor in neighbors:
            neighbor_pdb_id = self.atom_id_to_pdb_id[(neighbor.atomno(), neighbor.rsd())]
            adj_info[neighbor_pdb_id] = 1
        return adj_info
        
    # Prints the adjacency list
    def print_bond_adj_list(self):
        for atom_id in self.bond_adj_list:
            # Get atom id data
            bonded_atom_ids = self.bond_adj_list[atom_id]
            for bonded_atom_id in bonded_atom_ids:
                atom_res_index = atom_id.rsd()
                atom_index = atom_id.atomno()
                atom_residue = self.pose.residue(atom_res_index)
                atom_res_name = atom_residue.name()
                atom_name = atom_residue.atom_name(atom_index)
                # Get bonded neighbor's data
                bonded_atom_res_index = bonded_atom_id.rsd()
                bonded_atom_index = bonded_atom_id.atomno()
                bonded_atom_residue = self.pose.residue(bonded_atom_res_index)
                bonded_atom_res_name = bonded_atom_residue.name()
                bonded_atom_name = bonded_atom_residue.atom_name(bonded_atom_index)
                bonded_atom_pdb_index = self.atom_id_to_pdb_id[(bonded_atom_index, bonded_atom_res_index)]
                print("{} with Atom ID ({},{}) is bonded to {} with Atom ID({}, {}) and PDB_ID {} ".format(atom_name, atom_index, atom_res_index, bonded_atom_name, bonded_atom_index, bonded_atom_res_index, bonded_atom_pdb_index))
    
    # Debugging tool that checks to make sure each atom is in each other's bonded list
    def validate_bond_adj_list(self):
        valid = True
        count = 0
        # Loop through all atom ids
        for atom_id in self.bond_adj_list:
            # Loop through all atom neighbors
            bonded_atom_ids = self.bond_adj_list[atom_id]
            for bonded_atom_id in bonded_atom_ids:
                bonded_atomno = bonded_atom_id.atomno()
                bonded_res = bonded_atom_id.rsd()
                # Find the atom id address with the matching atomno and res
                bond_atom_id_key = None 
                for bonded_key in self.bond_adj_list:
                    if bonded_key.atomno() == bonded_atomno and bonded_key.rsd() == bonded_res:
                        bond_atom_id_key = bonded_key
                        break
                if bond_atom_id_key == None:
                    raise Exception("For some reason, a bonded atom is not in the bond adjacency list")
                atom_id_found = False
                # Check if atom id is in its neighbor's bonded list
                for bonded_bonded_atom_id in self.bond_adj_list[bond_atom_id_key]:
                    if bonded_bonded_atom_id.atomno() == atom_id.atomno() and bonded_bonded_atom_id.rsd() == atom_id.rsd():
                        atom_id_found = True
                        break
                # Print the atom that was not found in its neighbor's list
                if not atom_id_found:
                    count += 1
                    atom_res_index = atom_id.rsd()
                    atom_index = atom_id.atomno()
                    atom_residue = self.pose.residue(atom_res_index)
                    atom_res_name = atom_residue.name()
                    atom_name = atom_residue.atom_name(atom_index)
                    print(atom_id)
                    print(atom_name, atom_res_name, atom_index)
        return valid

    ################# SCORING FUNCTION #################
    
    # Retuns the energy of the current configuration
    # Set only_hbond to true if you only want the hbond energy
    def set_score_function(self, only_hbond=False):
        if only_hbond:
            sfxn = ScoreFunction()
            sfxn.set_weight(ScoreType.hbond_lr_bb, 1.0)
            sfxn.set_weight(ScoreType.hbond_sr_bb, 1.0)
            sfxn.set_weight(ScoreType.hbond_bb_sc, 1.0)
            sfxn.set_weight(ScoreType.hbond_sc, 1.0)
            return sfxn
        else:
            return get_score_function(True)
    
    # Gets the score of the current protein conformation (reward)
    def get_score(self):
        return self.score_function(self.pose)

    #################### GENERAL DIHEDRAL FUNCTIONS####################

    # Gets a list of TorsionIds
    # We can either get "all" angles, only "backbone" angles, or only "sidechain" angles  
    def get_torsion_ids(self, torsion_type="all"):
        torsion_ids = []
        for res_index in range(1, self.num_residues + 1):
            if torsion_type != "sidechain":
                torsion_ids.append(TorsionID(res_index, TorsionType.BB, 1))
                torsion_ids.append(TorsionID(res_index, TorsionType.BB, 2))
            if torsion_type != "backbone":
                for chi_index in range(1, self.pose.residue(res_index).nchi() + 1):
                    torsion_ids.append(TorsionID(res_index, TorsionType.CHI, chi_index))
        return np.array(torsion_ids)
    
    # Gets the torsion angles at a given torsion id
    def get_torsion_angles(self, torsion_ids):
        torsion_angles = np.zeros(len(torsion_ids))
        for index, torsion_id in enumerate(torsion_ids):
            res_index = torsion_id.rsd()
            torsion_type = torsion_id.type()
            torsion_num = torsion_id.torsion()
            #  Set the backbone angle
            if torsion_type == TorsionType.BB:
                if torsion_num == 1:
                    torsion_angles[index] = self.pose.phi(res_index)
                elif torsion_num == 2:
                    torsion_angles[index] = self.pose.psi(res_index)
                else:
                    raise Exception("Invalid torsion number for backbone")
            # Set the sidechain angle
            elif torsion_type == TorsionType.CHI:
                torsion_angles[index] = self.get_sidechain_angle(res_index, torsion_num)
            else:
                raise Exception("Invalid torsion type")
        return torsion_angles
    
    # Sets all the torsion ids in the list to the new torsions
    def set_torsion_ids(self, torsion_ids, new_torsions):
        # Loop through all the torsion ids
        for index, torsion_id in enumerate(torsion_ids):
            res_index = torsion_id.rsd()
            torsion_type = torsion_id.type()
            torsion_num = torsion_id.torsion()
            #  Set the backbone angle
            if torsion_type == TorsionType.BB:
                if torsion_num == 1:
                    self.pose.set_phi(res_index, new_torsions[index])
                elif torsion_num == 2:
                    self.pose.set_psi(res_index, new_torsions[index])
                else:
                    raise Exception("Invalid torsion number for backbone")
            # Set the sidechain angle
            elif torsion_type == TorsionType.CHI:
                self.pose.set_chi(torsion_num, res_index, new_torsions[index])
            else:
                raise Exception("Invalid torsion type")

    # Pertubs all torsions in list by torsion change                
    def perturb_torsion_ids(self, torsion_ids, torsion_changes):
        # Loop through all torsion ids
        for index, torsion_id in enumerate(torsion_ids):
            res_index = torsion_id.rsd()
            torsion_type = torsion_id.type()
            torsion_num = torsion_id.torsion()
            # Perturb backbone angles
            if torsion_type == TorsionType.BB:
                if torsion_num == 1:
                    self.pose.set_phi(res_index, self.pose.phi(res_index) + torsion_changes[index])
                elif torsion_num == 2:
                    self.pose.set_psi(res_index, self.pose.psi(res_index) + torsion_changes[index])
                else:
                    raise Exception("Invalid torsion number for backbone")
            # Perturb sidechain angle
            elif torsion_type == TorsionType.CHI:
                self.pose.set_chi(torsion_num, res_index, self.pose.chi(torsion_num, res_index) + torsion_changes[index])
            else:
                raise Exception("Invalid torsion type")
        return        

    ################# BACKBONE DIHEDRALS (PHI, PSI, OMEGA) #################

    # Returns an array of length 3 with the mainchain torsions
    # May need to change if we find out it does not trigger angle update on cart update
    def get_backbone_dihedrals(self, residue_index):
        return np.array(self.pose.residue(residue_index).mainchain_torsions())

    # Sets the phi angle of a residue (Uses 1-based indexing)
    def set_phi_angle(self, residue_index, new_angle):
        self.pose.set_phi(residue_index, new_angle)
        return

    # Sets the psi angle of a residue (Uses 1-based indexing)
    def set_psi_angle(self, residue_index, new_angle):
        self.pose.set_psi(residue_index, new_angle)
        return
    
    # Given a residue index, change the phi-psi angles of that residue
    # The dim of angle_change and cur_dihedrals are both vectors of length 2 
    # new_dihedrals = cur_dihedrals + angle_change 
    def perturb_backbone_dihedrals(self, residue_index, angle_change):
        residue = self.pose.residue(residue_index)
        #cur = [self.pose.phi(residue_index),  self.pose.psi(residue_index)]
        #new_phi_psi = np.append(cur + np.array(angle_change), self.pose.omega(residue_index))
        #residue.mainchain_torsions(Vector1(new_phi_psi))
        self.pose.set_phi(self.pose.phi(residue_index) + angle_change[0])
        self.pose.set_psi(self.pose.psi(residue_index) + angle_change[1])
        return
    
    # Sample backbone dihedrals from Dunbrack ramachanndran distribution
    # with bias if there is a proline before a residue 
    def sample_rama_backbone_dihedrals(self):
        RandomizeBBByRamaPrePro().apply(self.pose)
        return
    
    def sample_backbone_uniform(self):
        for res_index in range(1, self.num_residues + 1):
            new_phi = random.uniform(-180, 180)
            new_psi = random.uniform(-180, 180)
            self.set_phi_angle(res_index, new_phi)
            self.set_phi_angle(res_index, new_psi)
        return


    ################# SIDECHAIN DIHEDRALS (CHI) ROTAMERS #################

    # Gets the sidechain angle (i.e. chi) given the residue and angle index
    def get_sidechain_angle(self, residue_index, angle_index):
        return self.pose.chi(angle_index, residue_index)
    
    # Gets all the sidechain angles
    def get_sidechain_angles(self, res_index):
        return np.array(self.pose.residue(res_index).chi())

    # Sets one of the rotamer angles (i.e. side chain angle) of a specific residues
    # If one of the rotamers is not specified, we assume it remains the same
    def set_sidechain_angle(self, res_index, angle_index, new_angle):
        self.pose.set_chi(angle_index, res_index, new_angle)
        return
    
    # Perturbs all the sidechain angles of a residue by some angle change
    def perturb_sidechain_angles(self, res_index, angle_change):
        residue = self.pose.residue(res_index)
        # Store current chi angles
        cur_chis = residue.chi()
        # Only choose the residue angles that are needed
        res_angle_change = angle_change[:residue.nchi()]
        # Compute and set new chis
        new_chis = np.array(cur_chis) + res_angle_change
        self.set_sidechain_angles(res_index, new_chis)
        return
    
    # Sets all the sidechain angles in a residue
    def set_sidechain_angles(self, res_index, new_chis):
        np_new_chis = np.array(new_chis)
        residue = self.pose.residue(res_index)
        #assert residue.nchi() == len(new_chis)
        for chi_index in range(1, residue.nchi()+1):
            self.pose.set_chi(chi_index, res_index, np_new_chis[chi_index-1]) # the -1 is to account for 1-based indexing
        return
        
    # Setup for the rotamer set
    # 1. Discrete Backbone Dependent Rotamer Set
    # 2. Discrete Backbone Independent Rotamer Set
    # 3. Continuous Backbone Dependent Rotamer Set
    def setup_dep_rot_set(self, res_index, backbone_changed=False):
        rot_set = self.rot_set_fact.create_rotamer_set(self.pose)
        rot_set.set_resid(res_index)
        rot_set.build_rotamers(self.pose, self.score_function, self.pack_task, self.pack_graph)
        return rot_set
    
    # Sample backbone dependent rotamers (Takes a little longer than backbone independent)
    # Currently we assume that the backbone does not change
    def sample_bbdep_rotamers(self, backbone_changed=False):
        # Set each residue's rotamers
        for res_index in range(1, self.num_residues + 1):
            # Get the rotamer set for residue
            rot_set = self.setup_dep_rot_set(res_index)
            num_rotamers = rot_set.num_rotamers()
            # Sample random rotamer from rotamer set
            rand_rot_id = random.randint(1, num_rotamers)
            rotamers = rot_set.rotamer(rand_rot_id).chi()
            # Set the rotamer
            self.set_sidechain_angles(res_index, rotamers)
        return

    # Sets up the backbone independent rotamer set
    def setup_ind_rot_set(self):
        self.ind_rot_set = {}
        for res_index in range(1, self.num_residues + 1):
            residue_type = self.pose.residue_type(res_index)
            residue_type_name = str(residue_type.aa())
            # Skip if we already added this residue type
            if self.ind_rot_set.get(residue_type_name) != None:
                continue
            res_rotamers = rotamer_set.bb_independent_rotamers(residue_type)
            # Extract just the chi rotamers from the new residues
            rotamers = []
            for res_rot_index in range(1, len(res_rotamers) + 1):
                rotamers.append(res_rotamers[res_rot_index].chi())
            self.ind_rot_set[residue_type_name] = rotamers

    # Sample backbone indepedent rotamers    
    def sample_bbind_rotamers(self):
        for res_index in range(1, self.num_residues + 1):
            residue_type_name = str(self.pose.residue_type(res_index).aa())
            rot_set = self.ind_rot_set[residue_type_name]
            # Sample one of the rotamers
            rand_rot_id = random.randint(0, len(rot_set)-1)
            self.set_sidechain_angles(res_index, rot_set[rand_rot_id])
        return

    # Samples rotamers from a unifom distribution
    def sample_uniform_rotamers(self):
        for res_index in range(1, self.num_residues + 1):
            residue = self.pose.residue(res_index)
            num_rotamers = residue.nchi()
            # Skip if there are no rotamers
            if num_rotamers == 0:
                continue
            res_rotamers = []
            # Sample rotamers from a unifrom distribution [-180, 180]
            for i in range(num_rotamers):
                res_rotamers.append(random.uniform(-180, 180))
            self.set_sidechain_angles(res_index, res_rotamers)

    ################# HYDROGEN BONDS #################

    # Returns a dictionary mapping from an accceptor id to donor atom id 
    # used in the hydrogen bonding network
    def get_hbonds(self, rosetta=True):
        hbond_adj_list = {}
        # Generate the hbonds 
        self.hbond_set = self.pose.get_hbonds()
        # Loop through all hbonds 
        for hbond in self.hbond_set.hbonds():
            # Grab acceptor atom and residue index
            acc_res_index = hbond.acc_res()
            acc_atom_index = hbond.acc_atm()
            # Grab donor atom and residue index
            don_res_index = hbond.don_res()
            # For the donor atom, we want its parent not the actual hydrogen
            don_atom_index = self.pose.residue(don_res_index).atom_base(hbond.don_hatm())
            # Make atom id
            acc_atom_id = AtomID(acc_atom_index, acc_res_index) 
            don_atom_id = AtomID(don_atom_index, don_res_index) 
            hbond_adj_list[acc_atom_id] = don_atom_id
        return hbond_adj_list

    def print_hbond_adj_list(self, hbond_adj_list):
        for acc_atom_id in hbond_adj_list:
            don_atom_id = hbond_adj_list[acc_atom_id]
            # Grab acceptor atom info
            acc_res_index = acc_atom_id.rsd()
            acc_atom_index = acc_atom_id.atomno()
            acc_residue = self.pose.residue(acc_res_index)
            acc_res_name = acc_residue.name()
            acc_atom_name = acc_residue.atom_name(acc_atom_index)
            # Grab donor atom info
            don_res_index = don_atom_id.rsd()
            don_atom_index = don_atom_id.atomno()
            don_residue = self.pose.residue(don_res_index)
            don_res_name = don_residue.name()
            don_atom_name = don_residue.atom_name(don_atom_index)
            if(acc_residue.atom_type(acc_atom_index).is_virtual() or don_residue.atom_type(don_atom_index).is_virtual() ):
                print("VIRTUAL ATOM DETECTED")
            print("{} with Atom ID ({},{}) donated to a {} with Atom ID({}, {})".format(acc_atom_name, acc_atom_index, acc_res_index, don_atom_name, don_atom_index, don_res_index))

    # Prints the hbond network
    def print_hbond_network(self):
        if self.hbond_set == None:
            raise Exception("You must call get_hbonds to generate the hbond network before printing it")
        self.hbond_set.show(self.pose)
                    
    ################# ATOMWISE CHEMICAL FEATURES #################
    
    # Generates atom features using Rosetta
    # Each atom has the following 13 features in this order
    # 1. X Y Z (length 3)
    # 1. The atom element name (or one hot encoding of it (length 5))
    # 2. The atom the lj radius
    # 3. The atomic charge
    # 4. Boolean indicating whether atom could be donor
    # 5. Boolean indicating whether atom could be acceptor
    # 6. Boolean indicating whether atom is on the backbone
    def generate_chemical_features(self, one_hot=True, include_adj=True):
        chemical_features = []
        for atom_id in  self.atom_ids:
            chemical_feature = []
            # Grab cart coords
            xyz = self.pose.xyz(atom_id)
            xpos = xyz[0]
            ypos = xyz[1]
            zpos = xyz[2]
            # Grab atom index and residue
            atom_index = atom_id.atomno()
            residue = self.pose.residue(atom_id.rsd())
            atom_type = residue.atom_type(atom_index)
            # Grab all other chemical features
            lj_radius = atom_type.lj_radius()
            is_donor = atom_type.is_donor()
            is_acceptor = atom_type.is_acceptor()
            is_backbone = residue.atom_is_backbone(atom_index)
            atomic_charge = residue.atomic_charge(atom_index)
            #charge = residue.atomic_charge(atom_index)
            #hybrid = atom_type.hybridization()
            #atom_name = residue.atom_name(atom_index)
            
            # Add cart coords first
            chemical_feature = [xpos, ypos, zpos]
            # Add element 
            if one_hot:
                element_id = one_hot_map[atom_type.element()] # One hot encoding of element
                one_hot_element = [0] * NUM_ELEMENTS
                one_hot_element[element_id] = 1
                chemical_feature.extend(one_hot_element)
            else: 
                element = atom_type.element() # Element as string
                chemical_feature.append(element)
            # Add all other chemical features
            chemical_feature.extend([lj_radius, atomic_charge, is_donor, is_acceptor, is_backbone])
            #if self.incl_adj:
            #    bond_adj_info = self.get_bond_adj_for_atom(atom_id)
            #    #hbond_adj_info = self.get_hbond_adj_for_atom(atom_id)
            #    chemical_feature.extend(bond_adj_info)
            #    #chemical_feature.extend(hbond_adj_info)
            chemical_features.append(chemical_feature)
        # Convert chemical features to matrix
        chemical_features = np.array(chemical_features)
        return chemical_features
    
    # Resets the pose to original pose
    def reset_pose(self):
        self.pose = self.original_pose.clone()
    
    # Relaxes the protein by adjusting backbone dihedrals and sidechain dihedrals
    # to minimize the score function
    def relax_prot(self, max_iter=100):
        relax = FastRelax()
        relax.set_scorefxn(self.score_function)
        relax.max_iter(max_iter)
        relax.apply(self.pose)
        return
    
    # THIS DOES NOT WORK
    def pack_prot(self):
        # create a standard ScoreFunction
        scorefxn = get_fa_scorefxn() #  create_score_function_ws_patch('standard', 'score12')

        ############
        # PackerTask
        # a PackerTask encodes preferences and options for sidechain packing, an
        #    effective Rosetta methodology for changing sidechain conformations, and
        #    design (mutation)
        # a PackerTask stores information on a per-residue basis
        # each residue may be packed or designed
        # PackerTasks are handled slightly differently in PyRosetta
        ####pose_packer = PackerTask()    # this line will not work properly
        pose_packer = standard_packer_task(self.pose)
        # the pose argument tells the PackerTask how large it should be

        # sidechain packing "optimizes" a pose's sidechain conformations by cycling
        #    through (Dunbrack) rotamers (sets of chi angles) at a specific residue
        #    and selecting the rotamer which achieves the lowest score,
        #    enumerating all possibilities for all sidechains simultaneously is
        #    impractically expensive so the residues to be packed are individually
        #    optimized in a "random" order
        # packing options include:
        #    -"freezing" the residue, preventing it from changing conformation
        #    -including the original sidechain conformation when determining the
        #        lowest scoring conformation
        pose_packer.restrict_to_repacking()    # turns off design
        pose_packer.or_include_current(True)    # considers original conformation
        print( pose_packer )

        # packing and design can be performed by a PackRotamersMover, it requires
        #    a ScoreFunction, for optimizing the sidechains and a PackerTask,
        #    setting the packing and design options
        packmover = minimization_packing.PackRotamersMover(scorefxn, pose_packer)
        packmover.apply(self.pose)
        return