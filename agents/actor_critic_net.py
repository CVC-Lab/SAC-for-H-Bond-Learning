from custom_nn_modules.Graph_NNs import GraphConvolution, GraphAggregation, MLP
from utils import soft_update, hard_update, create_log_gaussian, logsumexp, batch_flat_to_tensors
import torch.nn.functional as F
from torch.distributions import Normal
import torch
import torch.nn as nn

LOG_SIG_MAX = 2
LOG_SIG_MIN = -20
epsilon = 1e-6

# Initialize Policy weights
def weights_init_(m):
    if isinstance(m, nn.Linear):
        nn.init.xavier_uniform_(m.weight, gain=1)
        nn.init.constant_(m.bias, 0)

class Actor(nn.Module):
    def __init__(self, conv_dim, node_dim, edge_dim, z_dim, action_space=None, dropout=0):
        super(Actor, self).__init__()
        # Store dims
        graph_conv_dim, aux_dim, linear_dim = conv_dim
        self.z_dim = z_dim
        
        self.gcn_layer = GraphConvolution(node_dim, graph_conv_dim, edge_dim, dropout)
        self.agg_layer = GraphAggregation(graph_conv_dim[-1], aux_dim, node_dim, dropout)
        self.mlp = MLP(aux_dim, linear_dim, nn.Tanh())
        self.last_layer = nn.Linear(linear_dim[-1], z_dim)

        # action rescaling
        if action_space is None:
            self.action_scale = torch.tensor(1.)
            self.action_bias = torch.tensor(0.)
        else:
            self.action_scale = torch.FloatTensor(
                (action_space.high - action_space.low) / 2.)
            self.action_bias = torch.FloatTensor(
                (action_space.high + action_space.low) / 2.)
        
        self.apply(weights_init_)
    
    def forward(self, state, tensor_shapes, hidden=None, activation=None):
        node, adj = batch_flat_to_tensors(state, tensor_shapes=tensor_shapes)
        adj = adj[:,:,:,1:].permute(0,3,1,2)
        input1 = torch.cat((hidden, node), -1) if hidden is not None else node 
        h = self.gcn_layer(input1, adj)
        input1 = torch.cat((h, hidden, node) if hidden is not None else (h, node), -1)
        h = self.agg_layer(input1, torch.tanh)
        h = self.mlp(h)
        h = self.last_layer(h)
        return h[:, :self.z_dim], h[:, self.z_dim:]
    
    def sample(self, state, tensor_shapes):
        mean, log_std = self.forward(state, tensor_shapes)
        std = log_std.exp()
        normal = Normal(mean, std)
        x_t = normal.rsample()  # for reparameterization trick (mean + std * N(0,1))
        y_t = torch.tanh(x_t)
        action = y_t * self.action_scale + self.action_bias
        log_prob = normal.log_prob(x_t)
        # Enforcing Action Bound
        log_prob -= torch.log(self.action_scale * (1 - y_t.pow(2)) + epsilon)
        log_prob = log_prob.sum(1, keepdim=True)
        mean = torch.tanh(mean) * self.action_scale + self.action_bias
        return action, log_prob, mean

    def to(self, device):
        self.action_scale = self.action_scale.to(device)
        self.action_bias = self.action_bias.to(device)
        return super(Actor, self).to(device)

# Each critic consists of two qnets
class Critic(nn.Module):
    def __init__(self, conv_dim, node_dim, edge_dim, z_dim, action_dim, num_nodes, input_action_dim, dropout=0):
        self.critic1 = QNet(conv_dim, node_dim, edge_dim, z_dim, action_dim, num_nodes, input_action_dim, dropout=0)
        self.critic2 = QNet(conv_dim, node_dim, edge_dim, z_dim, action_dim, num_nodes, input_action_dim, dropout=0)
        self.apply(weights_init_)
    
    def forward(self, state, action, tensor_shapes):
        x1 = self.critic1(state, action, tensor_shapes)
        x2 = self.critic1(state, action, tensor_shapes)
        return x1, x2

class QNet(nn.Module):

    # Structure: GraphConv Layer (GCL)-> Aggreation of previous GCL ->MLP
    # node_dim (int): dim of node feature
    # conv_dim ([int], int, int)):
    #   tuple containing hidden dims of each conv, output dim of aggregation, dim of last linear layer after GCN
    # edge_dim (int): dimension of edges
    # z_dim (int): Final layer for state processing
    # action_dim ([int]): linear_dims for MLPs processing action

    def __init__(self, conv_dim, node_dim, edge_dim, z_dim, action_dim, num_nodes, input_action_dim, dropout=0):
        super(QNet, self).__init__()
        graph_conv_dim, aux_dim, linear_dim = conv_dim
        self.node_dim = node_dim
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
    def forward(self, state, action, tensor_shapes, hidden=None, activation=None):
        node, adj = batch_flat_to_tensors(state, tensor_shapes=tensor_shapes)
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