from os import stat_result
import numpy as np
from agents.Base_Agent import Base_Agent
from mol_processors.RosettaProtein import RosettaProt
from utils import generate_one_hot_encoding, tensors_to_batch_flat, batch_flat_to_tensors
import torch
import torch.nn as nn
from utilities.OU_Noise import OU_Noise
from utilities.data_structures.Replay_Buffer import Replay_Buffer
from torch.optim import Adam
import torch.functional as F
from agents.actor_critic_agents.SAC import SAC
from custom_nn_modules.Graph_NNs import GraphConvolution, GraphAggregation, MLP

class Actor(nn.Module):
    def __init__(self, conv_dim, node_dim, edge_dim, z_dim, dropout=0):
        super(Actor, self).__init__()
        graph_conv_dim, aux_dim, linear_dim = conv_dim
        self.gcn_layer = GraphConvolution(node_dim, graph_conv_dim, edge_dim, dropout)
        self.agg_layer = GraphAggregation(graph_conv_dim[-1], aux_dim, node_dim, dropout)
        self.mlp = MLP(aux_dim, linear_dim, nn.Tanh())
        self.last_layer = nn.Linear(linear_dim[-1], z_dim)
    
    # adj is the adjacency matrix while node is the feature matrix
    # adj (batch, n, n, edge_dim)
    # tensor_shapes allows us unpack the state tensor to the right size
    def forward(self, state, tensor_shapes, hidden=None, activation=None):
        node, adj = batch_flat_to_tensors(state, tensor_shapes=tensor_shapes)
        adj = adj[:,:,:,1:].permute(0,3,1,2)
        input1 = torch.cat((hidden, node), -1) if hidden is not None else node 
        h = self.gcn_layer(input1, adj)
        input1 = torch.cat((h, hidden, node) if hidden is not None else (h, node), -1)
        h = self.agg_layer(input1, torch.tanh)
        h = self.mlp(h)
        
        return self.last_layer(h)

class Critic(nn.Module):

    # Structure: GraphConv Layer (GCL)-> Aggreation of previous GCL ->MLP
    # node_dim (int): dim of node feature
    # conv_dim ([int], int, int)):
    #   tuple containing hidden dims of each conv, output dim of aggregation, dim of last linear layer after GCN
    # edge_dim (int): dimension of edges
    # z_dim (int): Final layer for state processing
    # action_dim ([int]): linear_dims for MLPs processing action

    def __init__(self, conv_dim, node_dim, edge_dim, z_dim, action_dim, num_nodes, input_action_dim, dropout=0):
        super(Critic, self).__init__()
        graph_conv_dim, aux_dim, linear_dim = conv_dim
        self.total_node_dim = node_dim * num_nodes
        self.node_dim = node_dim
        self.num_nodes = num_nodes
        self.input_action_dim = input_action_dim
        # Process State
        self.gcn_layer = GraphConvolution(node_dim, graph_conv_dim, edge_dim, dropout)
        self.agg_layer = GraphAggregation(graph_conv_dim[-1], aux_dim, node_dim, dropout)
        self.mlp = MLP(aux_dim, linear_dim, nn.Tanh())
        self.last_state_layer = nn.Linear(linear_dim[-1], z_dim)
        # Processes action
        self.action_mlp = MLP(self.input_action_dim, action_dim, nn.Tanh())
        # Processes action and state
        self.last_layer = nn.Linear(action_dim[-1] + z_dim, 1)
    
    # adj is the adjacency matrix while node is the feature matrix
    def forward(self, state_action, tensor_shapes, hidden=None, activation=None):
        node, adj, action = batch_flat_to_tensors(state_action, tensor_shapes=tensor_shapes)
        # Process state        
        adj = adj[:,:,:,1:].permute(0,3,1,2)
        input1 = torch.cat((hidden, node), -1) if hidden is not None else node
        h = self.gcn_layer(input1, adj)
        input1 = torch.cat((h, hidden, node) if hidden is not None else (h, node), -1)
        h = self.agg_layer(input1, torch.tanh)
        h = self.mlp(h)
        h = self.last_state_layer(h)
        # Process action
        h2 = self.action_mlp(action)
        # Concatenate procesed state and action and process both
        h2 = torch.cat((h , h2))
        return self.last_layer(h2)
    
class SAC_Hbond(SAC):
    # Initializes the SAC hbond agent
    def __init__(self, config):
        Base_Agent.__init__(self, config)
        assert self.action_types == "CONTINUOUS", "Action types must be continuous. Use SAC Discrete instead for discrete actions"
        assert self.config.hyperparameters["Actor"]["final_layer_activation"] != "Softmax", "Final actor layer must not be softmax"
        # Extract hyperparams
        self.hyperparameters = config.hyperparameters["Actor_Critic_Agents"]
        self.node_dim = len(self.environment.features[0])
        self.num_nodes = len(self.environment.features)
        self.edge_dim = 1
        self.input_action_dim = len(self.environment.torsion_ids_to_change)
        # Initialize Critic
        crit_hyp = self.hyperparameters["Critic"]
        self.critic_local = Critic(crit_hyp["conv_dim"], self.node_dim, self.edge_dim, crit_hyp["z_dim"], crit_hyp["action_dim"], self.num_nodes, self.input_action_dim)
        self.critic_local_2 = Critic(crit_hyp["conv_dim"], self.node_dim, self.edge_dim, crit_hyp["z_dim"], crit_hyp["action_dim"], self.num_nodes, self.input_action_dim)
        self.critic_optimizer = torch.optim.Adam(self.critic_local.parameters(),
                                                 lr=self.hyperparameters["Critic"]["learning_rate"], eps=1e-4)
        self.critic_optimizer_2 = torch.optim.Adam(self.critic_local_2.parameters(),
                                                   lr=self.hyperparameters["Critic"]["learning_rate"], eps=1e-4)
        self.critic_target = Critic(crit_hyp["conv_dim"], self.node_dim, self.edge_dim, crit_hyp["z_dim"], crit_hyp["action_dim"], self.num_nodes, self.input_action_dim)
        self.critic_target_2 = Critic(crit_hyp["conv_dim"], self.node_dim, self.edge_dim, crit_hyp["z_dim"], crit_hyp["action_dim"], self.num_nodes, self.input_action_dim)
        # Copy params from local to target critic network
        Base_Agent.copy_model_over(self.critic_local, self.critic_target)
        Base_Agent.copy_model_over(self.critic_local_2, self.critic_target_2)
        # Initialize replay buffer
        self.memory = Replay_Buffer(self.hyperparameters["Critic"]["buffer_size"], self.hyperparameters["batch_size"],
                                    self.config.seed)
        # Initialize Action
        act_hyp = self.hyperparameters["Actor"]
        self.actor_local = Actor(act_hyp["conv_dim"], self.node_dim, self.edge_dim, act_hyp["z_dim"], self.num_nodes)
        self.actor_optimizer = torch.optim.Adam(self.actor_local.parameters(),
                                          lr=self.hyperparameters["Actor"]["learning_rate"], eps=1e-4)
        # Entropy Tuning
        self.automatic_entropy_tuning = self.hyperparameters["automatically_tune_entropy_hyperparameter"]
        if self.automatic_entropy_tuning:
            self.target_entropy = -torch.prod(torch.Tensor(self.environment.action_space.shape).to(self.device)).item() # heuristic value from the paper
            self.log_alpha = torch.zeros(1, requires_grad=True, device=self.device)
            self.alpha = self.log_alpha.exp()
            self.alpha_optim = Adam([self.log_alpha], lr=self.hyperparameters["Actor"]["learning_rate"], eps=1e-4)
        else:
            self.alpha = self.hyperparameters["entropy_term_weight"]

        # OU Noise
        self.add_extra_noise = self.hyperparameters["add_extra_noise"]
        if self.add_extra_noise:
            self.noise = OU_Noise(self.action_size, self.config.seed, self.hyperparameters["mu"],
                                  self.hyperparameters["theta"], self.hyperparameters["sigma"])

        self.do_evaluation_iterations = self.hyperparameters["do_evaluation_iterations"]
        return

    # Gets the title of the environment (i.e. Protein)
    def get_environment_title(self):
        return "SingleProtein"


    ########################## LOSS FUNCTIONS ###################
    def calculate_critic_losses(self, state_batch, action_batch, reward_batch, next_state_batch, mask_batch):
        """Calculates the losses for the two critics. This is the ordinary Q-learning loss except the additional entropy
         term is taken into account"""
        with torch.no_grad():
            next_state_action, next_state_log_pi, _ = self.produce_action_and_action_info(next_state_batch)
            qf1_next_target = self.critic_target(torch.cat((next_state_batch, next_state_action), 1))
            qf2_next_target = self.critic_target_2(torch.cat((next_state_batch, next_state_action), 1))
            min_qf_next_target = torch.min(qf1_next_target, qf2_next_target) - self.alpha * next_state_log_pi
            next_q_value = reward_batch + (1.0 - mask_batch) * self.hyperparameters["discount_rate"] * (min_qf_next_target)
        qf1 = self.critic_local(torch.cat((state_batch, action_batch), 1))
        qf2 = self.critic_local_2(torch.cat((state_batch, action_batch), 1))
        qf1_loss = F.mse_loss(qf1, next_q_value)
        qf2_loss = F.mse_loss(qf2, next_q_value)
        return qf1_loss, qf2_loss

    def calculate_actor_loss(self, state_batch):
        """Calculates the loss for the actor. This loss includes the additional entropy term"""
        action, log_pi, _ = self.produce_action_and_action_info(state_batch)
        qf1_pi = self.critic_local(torch.cat((state_batch, action), 1))
        qf2_pi = self.critic_local_2(torch.cat((state_batch, action), 1))
        min_qf_pi = torch.min(qf1_pi, qf2_pi)
        policy_loss = ((self.alpha * log_pi) - min_qf_pi).mean()
        return policy_loss, log_pi

    def calculate_entropy_tuning_loss(self, log_pi):
        """Calculates the loss for the entropy temperature parameter. This is only relevant if self.automatic_entropy_tuning
        is True."""
        alpha_loss = -(self.log_alpha * (log_pi + self.target_entropy).detach()).mean()
        return alpha_loss