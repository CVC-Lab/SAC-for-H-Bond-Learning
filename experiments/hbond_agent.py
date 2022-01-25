# Import system modules
import sys
import os
from datetime import datetime
from os.path import dirname, abspath
import itertools

from torch import onnx
from wandb import wandb_agent
sys.path.append(dirname(dirname(abspath(__file__))))
# Import Python Modules
import numpy as np
import wandb
import torch
# Import Pyrosetta
from pyrosetta import *
# Import custom modules
from environments.RosettaProtEnv import SingleProtEnv
from utilities.data_structures.replay_memory import ReplayMemory 
from agents.hbond_agent_new import SAC

#Inputs
# Define some constants
PDB_DIR = "../data/1BDD/"
OCTREE_DIR = "../Octree"
OUTPUT_DIR = "../data/1BDD/output"
# Define some parameters
pdb_id = "1bdd" # The pdb id of the input protein
dcd_file = ""
pdb_init_file = PDB_DIR + "1bdd.pdb"
pdb_out_file = PDB_DIR + "1cq0_pdb2pqr_charmm-outn.pdb"
psf_file = PDB_DIR + "1bdd_pdb2pqr_charmm.psf"

# Hyperparameters
hyperparameters =  {
    "output_model": "./results/hbond_agent_model",
    "device": torch.device("cuda" if torch.cuda.is_available() else "cpu"),
    "replay_size": 10000000,
    "seed": 123456,
    "start_steps": 1000, # Num steps until we start sampling
    "batch_size": 256,
    "updates_per_step": 1,
    "num_steps": 10000, # Max total number of time steps
    "eval": True,
    "discount_rate": 1, 
    "tau": 0.005,
    "alpha": 0.2,
    "lr": 0.0003,
    "policy": "Gaussian",
    "target_update_interval": 1,
    "env_name": "SingleProtEnv",
    "automatic_entropy_tuning": True,
    "output_pdb": "./results/AA/AA_eps_dynamite.pdb", # Output pdb of training episodes
    "output_pdb_test": "./results/AA/AA_eps_dynamite_test.pdb", # Output pdb of testing episodes
    "Env": {
        "torsions_to_change": "all", # all, backbone, or sidechain
        "adj_mat_type": "bonded", 
        "step_size": 0.5,
        "discount_rate": 1,
        "discount_rate_threshold": 100,
        "max_time_step": 300  # Maximum time step for each episode
    },
    "Actor": {
        "conv_dim":[[4 ,16], 128, [32, 64]], # Graph Convolution hidden dims, Aggregation output layer, MLP linear layers
    },
    "Critic": {
        "conv_dim":[[4 ,16], 128, [32, 64]],
        "z_dim": 8, # Output linear layer for processing state
        "action_dim": [32, 32] # MLP layer hidden dims to process actions
    }
}
# Extract hyperparams
output_pdb = hyperparameters['output_pdb'] 
replay_size = hyperparameters["replay_size"]
seed = hyperparameters["seed"]
start_steps = hyperparameters["start_steps"]
batch_size = hyperparameters["batch_size"]
updates_per_step = hyperparameters["updates_per_step"]
num_steps = hyperparameters["num_steps"]
policy = hyperparameters["policy"]
automatic_entropy_tuning = hyperparameters["automatic_entropy_tuning"]
env_name = hyperparameters["env_name"]
do_eval = hyperparameters["eval"]
output_model_file = hyperparameters["output_model"]
output_pdb_test = hyperparameters["output_pdb_test"]
device = hyperparameters["device"]

wandb_allowed = True
# Print Num Gpus
print("We have", torch.cuda.device_count(), "GPUs")
# Init Wandb
mode = "online" if wandb_allowed else "disabled"
project_name = "Sidechain Packing"
run_name = "run_" + datetime.now().strftime("%m:%d:%Y:%H:%M:%S")
group = "Dialanine"
notes = ""
wandb.login()
wandb.init(project=project_name, name=run_name, group=group, notes=notes, config=hyperparameters, mode=mode)
# Init Rosetta
init()
# Init environ
env = SingleProtEnv(hyperparameters["Env"], seq="AA")
# Set seeds
env.seed(seed)
env.action_space.seed(seed)
torch.manual_seed(seed)
np.random.seed(seed)
# Init Agent
agent = SAC(env.action_space, env, hyperparameters)
wandb.watch(agent.policy, log='gradients') # Log parameters and gradients of actor and critic
wandb.watch(agent.critic, log='gradients')
# Init Memory
memory = ReplayMemory(replay_size, seed)

# Training Loop
total_numsteps = 0
total_eval_numsteps = 0
eval_iepisode = 0
updates = 0
global_wandb_step = 0
# Store network topology in ONNX
dummy_state = torch.rand(1, agent.num_nodes * (agent.node_dim + agent.num_nodes)).to(device)
dummy_action = torch.rand(1, agent.input_action_dim).to(device)
agent.save_model('SingleProtEnv', suffix='.onnx', actor_input=(dummy_state), critic_input=(dummy_state, dummy_action))
wandb.save("./models/sac_actor_SingleProtEnv_.onnx")
wandb.save("./models/sac_critic_SingleProtEnv_.onnx")
del dummy_state
del dummy_action

# Loop through episodes
lowest_energy = np.finfo(np.float32).max
for i_episode in itertools.count(1):
    episode_reward = 0
    episode_steps = 0
    done = False
    # Reset environment
    new_output_pdb = output_pdb
    state = env.reset(new_output=new_output_pdb)
    init_energy = env.cur_score
    print("-------------------------------------")
    print("Init Energy: {}".format(init_energy))
    while not done:
        # Sample random action for a while
        if start_steps > total_numsteps:
            action = env.action_space.sample()  # Sample random action
        else:
            action = agent.select_action(state)  # Sample action from policy

        if len(memory) > batch_size:
            # Number of updates per step in environment
            for i in range(updates_per_step):
                # Update parameters of all the networks
                critic_1_loss, critic_2_loss, policy_loss, ent_loss, alpha = agent.update_parameters(memory, batch_size, updates)
                updates += 1
        
        # Transiton to next state 
        next_state, reward, done, _ = env.step(action) # Step

        # Ignore the "done" signal if it comes from hitting the time horizon.
        # (https://github.com/openai/spinningup/blob/master/spinup/algos/sac/sac.py)
        mask = 1 if episode_steps == env._max_episode_steps else float(not done)
        # Log Episode step info
        if total_numsteps % 25 == 0:
            print("-------------------------------------")
            print("Episode: {} | Episode Step: {} | Action: {} | Reward: {}, Energy: {}".format(i_episode, episode_steps, action, reward, env.prot.get_score()))
        memory.push(state, action, reward, next_state, mask) # Append transition to memory
        # Save Model
        if total_numsteps % 10 == 0:
            agent.save_model("SingleProtEnv")

        episode_steps += 1
        total_numsteps += 1
        episode_reward += reward
        state = next_state
    
    # Keep track of the lowest energy conformation found
    if env.cur_score < lowest_energy:
        wandb.run.summary["best_energy"] = env.cur_score
        wandb.run.summary["best_pose_step"] = total_numsteps 
        lowest_energy = env.cur_score
        
    print("Final Energy: {}".format(env.cur_score))
    wandb_log_dict = {"train_final_energy": env.cur_score, "train_init_energy": init_energy, "train_cum_reward": episode_reward,  "train_episode": i_episode, "total_numsteps": total_numsteps}
    # Log metrics
    wandb.log(wandb_log_dict)

    print("Episode: {}, total numsteps: {}, episode steps: {}, reward: {}".format(i_episode, total_numsteps, episode_steps, round(episode_reward, 2)))
    # Evaluate Agent and save to VMD
    if i_episode % 10 == 0 and do_eval is True:
        avg_reward = 0.
        avg_energy_change = 0.
        episodes = 10
        for _  in range(episodes):
            state = env.reset(new_output=output_pdb_test)
            init_energy = env.cur_score
            episode_reward = 0
            done = False
            while not done:
                action = agent.select_action(state, evaluate=True)

                next_state, reward, done, _ = env.step(action)
                episode_reward += reward

                state = next_state
            avg_reward += episode_reward
            avg_energy_change += env.cur_score - init_energy
            eval_iepisode += 1
            # Log per episode metrics
            wandb.log({"eval_episode": eval_iepisode, "eval_init_energy": init_energy, "eval_final_energy": env.cur_score, "eval_cum_reward": episode_reward})
        avg_reward /= episodes
        avg_energy_change /= episodes


        print("----------------------------------------")
        print("Test Episodes: {}, Avg. Reward: {}".format(episodes, round(avg_reward, 2)))
        print("----------------------------------------")
        # Log Average test metrics
        wandb.log({"eval_delta_energy": avg_energy_change, "eval_avg_reward": avg_reward})
        
    # Stop when we reach the max total steps
    if total_numsteps > num_steps:
        break
    
    # Save model one last time before finished
    agent.save_model("SingleProtEnv")
    global_wandb_step += 1

env.close()
