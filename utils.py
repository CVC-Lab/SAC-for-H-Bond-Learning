import numpy as np
import torch
import os
from datetime import datetime as dt
import shutil
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
import math
import MDAnalysis as mda
import torch

def create_log_gaussian(mean, log_std, t):
    quadratic = -((0.5 * (t - mean) / (log_std.exp())).pow(2))
    l = mean.shape
    log_z = log_std
    z = l[-1] * math.log(2 * math.pi)
    log_p = quadratic.sum(dim=-1) - log_z.sum(dim=-1) - 0.5 * z
    return log_p

def logsumexp(inputs, dim=None, keepdim=False):
    if dim is None:
        inputs = inputs.view(-1)
        dim = 0
    s, _ = torch.max(inputs, dim=dim, keepdim=True)
    outputs = s + (inputs - s).exp().sum(dim=dim, keepdim=True).log()
    if not keepdim:
        outputs = outputs.squeeze(dim)
    return outputs

def soft_update(target, source, tau):
    for target_param, param in zip(target.parameters(), source.parameters()):
        target_param.data.copy_(target_param.data * (1.0 - tau) + param.data * tau)

def hard_update(target, source):
    for target_param, param in zip(target.parameters(), source.parameters()):
        target_param.data.copy_(param.data)

# Normalizes an angle in between -180 and 180 degrees or -pi and pi radians
# in_degrees(bool): if True, angle is in degrees
def normalize_angle(angle, in_degrees=True):
    if in_degrees:
        angle = angle % 360
        angle = (angle + 360) % 360
        if (angle > 180):
            angle -= 360
    else:
        angle = angle % 2*math.pi
        angle = (angle + 2*math.pi) % (2*math.pi)
        if (angle > math.pi):
            angle -= 2*math.pi
    
# Converts a list of torch tensors to one long unbatched flat tensor
# tensors ([tensors]): iterable of tensors 
def tensors_to_flat(tensors):
    flat_tensors = []
    # Flatten each tensor
    for tensor in tensors:
        flat_tensors.append(torch.flatten(tensor))
    # Concatenate the flattened tensors
    return torch.cat(flat_tensors) 

# Converts one long flat batched tensor to multiple reshaped unbatched tensors
# The product of all the tensor shapes must match the length of flat tensor
# flat_tensor (tensor): flat batched tensor
# tensor_shapes ([[int]]): list of tensor shapes
def flat_to_tensors(flat_tensor, tensor_shapes):
    start_dim = 0
    unflattened_tensors = [] 
    # Loop through all shapes
    for tensor_shape in tensor_shapes:
        total_dim = np.prod(np.array(tensor_shape)) # Find product of dims 
        end_dim = start_dim + total_dim 
        tensor = flat_tensor[start_dim:end_dim] # Grabs correct part of the flat tensor
        unflattened_tensors.append(torch.reshape(tensor, tensor_shape))
        start_dim += total_dim
    return tuple(unflattened_tensors)

# Converts a list of torch tensors to one batched long flat tensor 
# tensors ([tensors]): iterable of tensors
def tensors_to_batch_flat(tensors):  
    flat_tensors = []
    # Flatten each tensor
    for tensor in tensors:
        flat_tensors.append(torch.flatten(tensor, start_dim=1))
    # Concatenate the flattened tensors
    return torch.cat(flat_tensors, 1) 
 
# Converts one long flat batched tensor to multiple reshaped batched tensors
# The product of all the tensor shapes must match the length of flat tensor
# flat_tensor (tensor): flat batched tensor
# tensor_shapes ([[int]]): list of tensor shapes
def batch_flat_to_tensors(flat_tensor, tensor_shapes):
    start_dim = 0
    unflattened_tensors = [] 
    # Loop through all shapes
    for tensor_shape in tensor_shapes:
        total_dim = np.prod(np.array(tensor_shape)) # Find product of dims 
        end_dim = start_dim + total_dim 
        tensor = flat_tensor[:, start_dim:end_dim] # Grabs correct part of the flat tensor
        unflattened_tensors.append(torch.reshape(tensor, (-1,) + tensor_shape))
        start_dim += total_dim
    return tuple(unflattened_tensors)

# Appends a new pdb file to current pdb file
# new_pdb_file (str): path to new pdb file
# cur_pdb_file (str): path to current pdb file
def append_pdb(new_pdb, cur_pdb_file):
    with open(cur_pdb_file, "a") as writer:
        new_pdb_file = open(new_pdb, "r")
        for line in new_pdb_file:
            # Only write the ATOM coordinate lines and END
            if line.startswith("ATOM") or line.startswith("TER") or line.startswith("END"):
                writer.write(line)
        writer.write("END\n")
        new_pdb_file.close()
        
# Combine multiple PDBsinto one pdb file
def combine_pdbs(pdb_files, output_file):
    with open(output_file, "w") as writer:
        for pdb_file_name in pdb_files:
            # Read in pdb file
            pdb_file = open(pdb_file_name, "r")
            # Read line by line
            for line in pdb_file:
                # Only write the ATOM coordinate lines and END
                if line.startswith("ATOM") or line.startswith("TER") or line.startswith("END"):
                    writer.write(line)
            pdb_file.close()
    return

# Converts "arr" to numpy array and saves it as text or NPY format
def write_array(arr, output_file=None, file_type=None):
    np_array = np.array(arr).round(decimals=3)
    if output_file == None:
        output_file = "test.txt"
    # Save as binary
    if file_type == "binary":
        np.save(output_file, np_array)
    # Save as text
    elif file_type == None or file_type == 'text':
        np.savetxt(output_file, np_array, fmt="%.3f")
    else:
        raise Exception("file_type must be either 'binary' or 'text'")
    return

# Creates a one hot encoding from an array of objects
# arr (iterable): iterable of objects
def generate_one_hot_encoding(arr):
    # Determine number of unique objects
    np_arr = np.array(arr)
    unique_objects = np.unique(np_arr)
    # Create mapping between object and one hot encoding 
    one_hot_mapping = {}
    for index, unique_obj in enumerate(unique_objects):
        one_hot_mapping[unique_obj] = index
    # Create one enoding of arr
    one_hot_encoding = np.zeros((len(arr), len(unique_objects)))
    for index, obj in enumerate(arr):
        one_hot_index = one_hot_mapping[obj]
        one_hot_encoding[index, one_hot_index] = 1  
    return one_hot_encoding

# Converts chain sequences in a dictionary into one sequence for MODELLER
# Params:
# seq_dict (dict): dictionary where chain ids are keys and sequences are values
# header (str): the user specified header of FASTA sequence
# pretty (bool): flag to place line breaks every 75 characters
def sum_seq_in_dict (seq_dict, header="", pretty=True):
    result = ""
    # Add up all the chain sequences
    for key in seq_dict:
        result += seq_dict[key]
    # Remove trailing "\" char
    result = result[:-1] + "*"
    # Insert line breaks if pretty
    if pretty:
        result = insert_linebreaks (result)
    return header + result

# Inserts line breaks at specified intervals
# Params:
# string (string): the input string
# interval (int): the interval to insert line breaks 
def insert_linebreaks(string, interval=75):
    lines = []
    for i in range (0, len(string), interval):
        lines.append (string[i:i+interval])
    return '\n'.join(lines)

# Replaces a char in a string at given index
def str_replace (old_string, new_string, index):
    if index <= 0:
        return new_string + old_string
    if index >= len(old_string):
        return old_string + new_string
    
    return old_string[:index] + new_string + old_string[index + 1:]

# Insert a char
def str_insert (old_string, new_string, index):
    if index <= 0:
        return new_string + old_string
    if index >= len(old_string):
        return old_string + new_string
    
    return old_string[:index] + new_string + old_string[index:]

# Prints the parameters of a model
def print_params(model):
    i = 0
    print("Parameter List:")
    for name, param in model.named_parameters():
            trained = "NOT be trained"
            if param.requires_grad:
                trained = "be trained"
            print ("Param ", name, " will", trained)
            i += 1
    print ("There are ", i, "parameters")

# Implements gaussian cross entropy
def gauss_cross_entropy(mu1, var1, mu2, var2):
    term0 = 1.8378770664093453
    term1 = torch.log(var2)
    term2 = (var1 + mu1 ** 2 - 2 * mu1 * mu2 + mu2 ** 2) / var2

    cross_entropy = -0.5 * (term0 + term1 + term2)

    return cross_entropy

# Checks if the tensor has a nan value
def check_nan (tensor):
    if(torch.sum (torch.isnan (tensor)) > 0):
        return True
    return False 

if __name__ == '__main__':
    print ("No compile errors")