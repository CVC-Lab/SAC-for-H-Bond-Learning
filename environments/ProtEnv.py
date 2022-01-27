import numpy as np
import gym
from gym import spaces
from mol_processors.Protein import Prot


class ProtEnv(gym.Env):

    def __init__(self, hyperparameters, pdb_path: str, top_file: str):
        # Load hyperparameters
        self.min_torsion_change = hyperparameters["min_torsion_angle"]
        self.max_torsion_change = hyperparameters["max_torsion_angle"]
        self.active_torsion_type = hyperparameters["active_torsion_type"]
        # Create protein class
        self.protein = Prot(pdb_path, top_file)
        self.torsion_indices = self.protein.get_torsion_indices(self.active_torsion_type)
        # Set active 

        # Action space for protein are the torsion angles
        self.action_space = spaces.Box(
            low=self.min_torsion_change,
            high=self.max_torsion_change,
            shape=(len(self.torsion_indices),),
            dtype=np.float32
        )
        # Observation space for protein
        # TODO: 
        self.observation_space = spaces.Box(
            low=self.low_state,
            high=self.high_state,
            shape=(num_atoms, num_features + num_atoms),
            dtype=np.float32
        )
        return
    
    def bond_length_energy(self):
        """
        Computes a harmonic bond length energy term
        """
        return 0
    
    def bond_angle_energy(self):
        """
        Computes a harmonic bond angle energy term
        """
        return 0
    
    def torsion_angle_energy(self):
        """
        Computes a torsion angle energy term
        """
        return 0
    
    def nonbonded_energy(self):
        """
        Computes a nonbonded energy term (i.e. Coulombic + Lennard-Jones)
        """
        return 0

    def step(self, action):
        state = None
        reward = None
        done = False
        return state, reward, done, {}
    
    def reset(self):
        state = None
        return state
    
    def get_reward(self):
        return
    