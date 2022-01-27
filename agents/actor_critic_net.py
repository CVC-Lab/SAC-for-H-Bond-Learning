from utils import soft_update, hard_update, create_log_gaussian, logsumexp, batch_flat_to_tensors
import torch.nn.functional as F
from torch.distributions import Normal
import torch
import torch.nn as nn
from custom_nn_modules import MLP


LOG_SIG_MAX = 2
LOG_SIG_MIN = -20
epsilon = 1e-6

# Initialize Policy weights
def weights_init_(m):
    if isinstance(m, nn.Linear):
        nn.init.xavier_uniform_(m.weight, gain=1)
        nn.init.constant_(m.bias, 0)

class Actor(nn.Module):
    def __init__(self, state_dim, action_dim, action_space=None, dropout=0):
        super(Actor, self).__init__()
        self.state_dim = state_dim
        self.num_actions = action_dim
        layers = [4, 4, 6, 8]

        self.net = MLP(state_dim, layers, activation=nn.Relu(), dropout=0)
        self.last_layer = nn.Linear(layers[-1], self.num_actions * 2)

        # Action scaling factor and bias (default none)
        if action_space is None:
            self.action_scale = torch.tensor(1.)
            self.action_bias = torch.tensor(0.)
        else:
            self.action_scale = torch.FloatTensor(
                (action_space.high - action_space.low) / 2.)
            self.action_bias = torch.FloatTensor(
                (action_space.high + action_space.low) / 2.)
        
        self.apply(weights_init_)
    
    def forward(self, state, hidden=None, activation=None):
        # Apply Linear layers
        h = self.net(state)
        h = self.last_layer(h)
        # Partition the output to mean and log std
        action_size = int(self.num_actions*2)
        mu = h[:, :action_size]
        log_std =  h[:, action_size:]
        log_std = torch.clamp(log_std, min=LOG_SIG_MIN, max=LOG_SIG_MAX)
        return mu, log_std
    
    def sample(self, state):
        mean, log_std = self.forward(state)
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
    def __init__(self, state_dim, action_dim, dropout=0):
        super(Critic, self).__init__()
        self.critic1 = QNet(state_dim, action_dim, dropout=0)
        self.critic2 = QNet(state_dim, action_dim, dropout=0)
        self.apply(weights_init_)
    
    def forward(self, state, action):
        x1 = self.critic1(state, action, self.tensor_shapes)
        x2 = self.critic1(state, action, self.tensor_shapes)
        return x1, x2

class QNet(nn.Module):

    def __init__(self, state_dim, action_dim, dropout=0):
        super(QNet, self).__init__()
        self.state_dim = state_dim
        self.num_actions = action_dim
        state_net_layers = [4, 4, 4]
        final_state_dim = 2
        act_net_layers = [4, 4, 4]
        final_action_dim = 2

        # Process state
        self.state_net = MLP(state_dim, state_net_layers, activation=nn.Relu(), dropout=0)
        self.last_state_layer = nn.Linear(state_net_layers[-1], final_state_dim)

        # Processes action
        self.action_net = MLP(action_dim, act_net_layers, activation=nn.Relu(), dropout=0)
        self.last_action_layer = nn.Linear(act_net_layers[-1], final_action_dim)

        # Process both action and state to only scalar representing the action-value
        self.last_layer = nn.Linear(action_dim[-1] + final_state_dim, 1)
    
    # adj is the adjacency matrix while node is the feature matrix
    def forward(self, state, action, tensor_shapes, hidden=None, activation=None):
        # Process state
        h = self.state_net(state)
        h = self.last_state_layer(h)
        # Process action
        h2 = self.action_net(action)
        h2 = self.last_action_layer(h2)
        # Concatenate processed state and action and compute Q-value
        h2 = torch.cat((h , h2), 1)
        return self.last_layer(h2)