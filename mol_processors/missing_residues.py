import os
from collections import defaultdict

from modeller import *
from modeller.automodel import *
from utils import  sum_seq_in_dict, str_insert
from mol_processors.pdb_processor import download_pdb, parse_pdb_header, res_3to1
from Bio.PDB import *


# TODO: Test missing residues module

# Fills in missing residues of pdb file using MODELLER
# The files will be ouputted to your current directory
# If the file path of the pdb file is not specified, a new pdb is downloaded in the current directory
# Params:
# pdb_id (str): The pdb id of the PDB file
# file_path (str): File path of the PDB file (NOTE THAT IT MUST START WITH "pdb_id")
# read_het (bool): Flag to also read in the HETATMS 
def fill_missing_residues(pdb_id, file_path=None,read_het=False):
    # Download the pdb it is not downloaded
    if file_path == None:
        download_pdb (pdb_id)
        file_path = os.getcwd() + "/" + pdb_id.lower() + ".pdb"
    # Check first if there are any missing residues
    head = parse_pdb_header (file_path)
    if not head['has_missing_residues']:
        print ("There are not missing residues")
        return
    # Generate full and missing sequence dictionaries
    full_seq, miss_seq = create_full_and_missing_seq (pdb_id, file_path, head=head, read_het=read_het)
    if (miss_seq == None):
        print ("There are too many consecutive missing residues")
        return
    # Set up MODELLER environment to read heteroatoms
    #log.verbose()
    env = environ()
    env.io.hetatm = read_het
    env.io.atom_files_directory = [os.path.dirname(file_path)]
    # Create a model of pdb id
    m = model(env, file='6gyp')
    aln = alignment(env), create_full_and_missing_seq
    with open (seq_file + ".seq", "r") as f:
        i = 0
        for line in f:
            if i <= 2:
                header.append(line)
            else:
                break
            i += 1
    # Write the align file containing the full sequence and the summed sequences
    align_path = os.path.dirname(file_path) + "/" + pdb_id + "_align.ali"
    f = open (align_path, "w")   
    f.write(sum_seq_in_dict (miss_seq, header=header[0] + header[1] + header[2], pretty=True))
    f.write("\n")
    f.write(sum_seq_in_dict (full_seq, header=header[0] + header[1][:-1] + "_fill\n" + "sequence:::::::::\n",pretty=True))
    f.close()

    # Perform the loop modeling
    a = loopmodel(env, alnfile = align_path,
              knowns = pdb_id, sequence = pdb_id + "_fill")
    a.starting_model= 1
    a.ending_model  = 1
    a.loop.starting_model = 1
    a.loop.ending_model   = 2
    a.loop.md_level       = refine.fast

    a.make()
    return


# Creates two dictionaries where each key refers to the chain id and each
# value refers to its amino acid sequence. In one dictionary, the full 
# sequence of each chain will be stored and in the other each missing residue
# is represented as a "-" char
# Params:
# pdb_id (str): The pdb id of the PDB file
# file_path (str): File path of the PDB file (NOTE THAT IT MUST START WITH "pdb_id")
# read_het (bool): Flag to also read in the HETATMS  
def create_full_and_missing_seq (pdb_id, file_path, head=None, read_het=False):
    
    # Extract header information if not given
    if head is None:
        head = parse_pdb_header (file_path)
    # Extract chains
    structure = PDBParser().get_structure(pdb_id.upper(), file_path)
    chains = [each.id for each in structure.get_chains()]
    print (chains)
    # Make a dictionary for each chain
    missing_res = defaultdict(list)
    # Fill each chain's list with missing residue indices and missing residue names
    for residue in head['missing_residues']:
        if residue["chain"] in chains:
            seq_num = residue["ssseq"]
            res_code = res_3to1[residue["res_name"]]
            missing_res[residue["chain"]].append({"seq_num": seq_num, "res_code": res_code})
        
    # Collect the full and missing sequences
    full_seq = {}
    miss_seq = {}
    miss_index = 0
    # Grab the original pdb sequence by chain
    for chain in structure[0]:
        orig_seq = ""
    # Extract each chain from the pdb file 
        print (">Chain", chain.id, "\n", ''.join(orig_seq))
        print ()
        # Try to get missing res list for chain
        miss_indices = missing_res.get (chain.id)
        # Case for no missing residues in chain
        if (miss_indices == None):
            print ("No missing residues for chain ", chain)
            full_seq[chain.id] = orig_seq + "/"
            miss_seq[chain.id] = miss_seq + "/"
        # Case for missing residues
        else:
            print ("There are", len(miss_indices), "missing residues in chain", chain.id)
            # Replace all 
            temp_full = orig_seq
            temp_miss = orig_seq
            # Add in missing residues or gap characters
            for res in miss_indices:
                seq_index = res["seq_num"] - 1
                temp_full = str_insert (temp_full, res["res_code"], seq_index)
                temp_miss = str_insert (temp_miss, "-", seq_index)
            # Add the chain sequence
            full_seq[chain.id] = temp_full + "/"
            miss_seq[chain.id] = temp_miss + "/"
    # Extract each chain from the pdb file 
    if not below_max_missing_res (miss_seq):
        return None, None
    return full_seq, miss_seq

# Checks if a sequence does not have a consecutive missing residues larger than
# a threshold
# Params:
# seq (dict): the dictonary of amino acid sequence separated by chain
# threshold (int): the max number of residues allowed
def below_max_missing_res (seq, threshold=100):
    for chain in seq:
        num_consec = 0
        for char in seq[chain]:
            if char == "-":
                num_consec += 1
                if num_consec > threshold:
                    print ("Chain ", chain, "had ", num_consec, "consecutive missing residues")
                    return False
            else:
                num_consec = 0
    return True

# Retuirns true if a REMARK line in PDB says there are missing residues
# pdb (string): path to pdb file
def missing_residue(pdb):
    for line in open(pdb):
        list = line.split()
        id = list[0]
        if id == 'REMARK':
            if list[1] == '465':
                return True
    return False

