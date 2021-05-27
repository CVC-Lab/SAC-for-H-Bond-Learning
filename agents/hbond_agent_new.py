import os
import torch
import torch.nn.functional as F
from torch.optim import Adam
from utils import soft_update, hard_update, create_log_gaussian, logsumexp
from custom_nn_modules.Graph_NNs import GraphConvolution, GraphAggregation, MLP
from agents.replay_memory import ReplayMemory
from agents.actor_critic_net import Actor, Critic

class SAC(object):
    def __init__(self, num_inputs, action_space, hyperparams):

        # Hyperparametres
        self.gamma = hyperparams["gamma"]
        self.tau = hyperparams["tau"]
        self.alpha = hyperparams["alpha"]
        self.lr = hyperparams["lr"]
        self.policy_type = hyperparams["policy_type"]
        self.target_update_interval = hyperparams["target_update_interval"]
        self.automatic_entropy_tuning = hyperparams["automatic_entropy_tuning"]
        self.environment = hyperparams["env"]
        # Layer Hyperparams
        self.node_dim = len(self.environment.features[0])
        self.num_nodes = len(self.environment.features)
        self.edge_dim = 1
        self.input_action_dim = len(self.environment.torsion_ids_to_change)
        crit_hyp = hyperparams["Critic"]
        act_hyp = hyperparams["Actor"]
        self.tensor_shapes = [(self.node_dim, self.num_nodes), (self.num_nodes, self.num_nodes)]
        
        self.device = torch.device("cuda" if hyperparams["cuda"] else "cpu")

        self.critic = Critic(crit_hyp["conv_dim"], self.node_dim, self.edge_dim, crit_hyp["z_dim"], crit_hyp["action_dim"], self.num_nodes, self.input_action_dim)
        self.critic_optim = Adam(self.critic.parameters(), lr=self.lr)

        self.critic_target = Critic(crit_hyp["conv_dim"], self.node_dim, self.edge_dim, crit_hyp["z_dim"], crit_hyp["action_dim"], self.num_nodes, self.input_action_dim)
        hard_update(self.critic_target, self.critic)

        if self.policy_type == "Gaussian":
            # Target Entropy = −dim(A) (e.g. , -6 for HalfCheetah-v2) as given in the paper
            if self.automatic_entropy_tuning is True:
                self.target_entropy = -torch.prod(torch.Tensor(action_space.shape).to(self.device)).item()
                self.log_alpha = torch.zeros(1, requires_grad=True, device=self.device)
                self.alpha_optim = Adam([self.log_alpha], lr=self.lr)

            self.policy = Actor(act_hyp["conv_dim"], self.node_dim, self.edge_dim, act_hyp["z_dim"], self.num_nodes)
            self.policy_optim = Adam(self.policy.parameters(), lr=self.lr)

        #else:
        #    self.alpha = 0
        #    self.automatic_entropy_tuning = False
        #    self.policy = DeterministicPolicy(num_inputs, action_space.shape[0], args.hidden_size, action_space).to(self.device)
        #    self.policy_optim = Adam(self.policy.parameters(), lr=args.lr)


    # Selects an action
    def select_action(self, state, evaluate=False):
        state = torch.FloatTensor(state).to(self.device).unsqueeze(0)
        if evaluate is False:
            action, _, _ = self.policy.sample(state, self.tensor_shapes)
        else:
            _, _, action = self.policy.sample(state, self.tensor_shapes)
        return action.detach().cpu().numpy()[0]

    # Updates the parameters of the actor and critic
    def update_parameters(self, memory, batch_size, updates):
        # Sample a batch from memory
        state_batch, action_batch, reward_batch, next_state_batch, mask_batch = memory.sample(batch_size=batch_size)

        state_batch = torch.FloatTensor(state_batch).to(self.device)
        next_state_batch = torch.FloatTensor(next_state_batch).to(self.device)
        action_batch = torch.FloatTensor(action_batch).to(self.device)
        reward_batch = torch.FloatTensor(reward_batch).to(self.device).unsqueeze(1)
        mask_batch = torch.FloatTensor(mask_batch).to(self.device).unsqueeze(1)

        with torch.no_grad():
            next_state_action, next_state_log_pi, _ = self.policy.sample(next_state_batch, self.tensor_shapes)
            qf1_next_target, qf2_next_target = self.critic_target(next_state_batch, next_state_action)
            min_qf_next_target = torch.min(qf1_next_target, qf2_next_target) - self.alpha * next_state_log_pi
            next_q_value = reward_batch + mask_batch * self.gamma * (min_qf_next_target)
        qf1, qf2 = self.critic(state_batch, action_batch)  # Two Q-functions to mitigate positive bias in the policy improvement step
        qf1_loss = F.mse_loss(qf1, next_q_value)  # JQ = 𝔼(st,at)~D[0.5(Q1(st,at) - r(st,at) - γ(𝔼st+1~p[V(st+1)]))^2]
        qf2_loss = F.mse_loss(qf2, next_q_value)  # JQ = 𝔼(st,at)~D[0.5(Q1(st,at) - r(st,at) - γ(𝔼st+1~p[V(st+1)]))^2]
        qf_loss = qf1_loss + qf2_loss

        self.critic_optim.zero_grad()
        qf_loss.backward()
        self.critic_optim.step()

        pi, log_pi, _ = self.policy.sample(state_batch, self.tensor_shapes)

        qf1_pi, qf2_pi = self.critic(state_batch, pi)
        min_qf_pi = torch.min(qf1_pi, qf2_pi)

        policy_loss = ((self.alpha * log_pi) - min_qf_pi).mean() # Jπ = 𝔼st∼D,εt∼N[α * logπ(f(εt;st)|st) − Q(st,f(εt;st))]

        self.policy_optim.zero_grad()
        policy_loss.backward()
        self.policy_optim.step()

        if self.automatic_entropy_tuning:
            alpha_loss = -(self.log_alpha * (log_pi + self.target_entropy).detach()).mean()

            self.alpha_optim.zero_grad()
            alpha_loss.backward()
            self.alpha_optim.step()

            self.alpha = self.log_alpha.exp()
            alpha_tlogs = self.alpha.clone() # For TensorboardX logs
        else:
            alpha_loss = torch.tensor(0.).to(self.device)
            alpha_tlogs = torch.tensor(self.alpha) # For TensorboardX logs


        if updates % self.target_update_interval == 0:
            soft_update(self.critic_target, self.critic, self.tau)

        return qf1_loss.item(), qf2_loss.item(), policy_loss.item(), alpha_loss.item(), alpha_tlogs.item()

    # Save model parameters
    def save_model(self, env_name, suffix="", actor_path=None, critic_path=None):
        if not os.path.exists('models/'):
            os.makedirs('models/')

        if actor_path is None:
            actor_path = "models/sac_actor_{}_{}".format(env_name, suffix)
        if critic_path is None:
            critic_path = "models/sac_critic_{}_{}".format(env_name, suffix)
        print('Saving models to {} and {}'.format(actor_path, critic_path))
        torch.save(self.policy.state_dict(), actor_path)
        torch.save(self.critic.state_dict(), critic_path)

    # Load model parameters
    def load_model(self, actor_path, critic_path):
        print('Loading models from {} and {}'.format(actor_path, critic_path))
        if actor_path is not None:
            self.policy.load_state_dict(torch.load(actor_path))
        if critic_path is not None:
            self.critic.load_state_dict(torch.load(critic_path))


