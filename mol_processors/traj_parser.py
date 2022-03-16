import MDAnalysis as mda
import numpy as np

# Grabs coordinates from a file and returns a numpy array of coords.
# Optionally save coords to pdb or npy file
# 
# coord_path (str): path to the trajectory or pdb file
# top_path (str): path to the topology file
# file_type (str): type of coord file (dcd or pdb)
# save_pdbs (bool): if True, saves the coordinates to a pdb file
# Return: N x 3 Numpy array of coords
def get_coords(coord_path, top_path, file_type="dcd", save_pdbs=False, save_np=False, np_file="prot_coords.npy"):
    u = mda.Universe(top_path, coord_path)
    protein = u.select_atoms("protein")
    result = []
    
    # Case if the file is a trajectory file
    if file_type == "dcd":
        # Add the postions of the protein at each timestep to result
        for ts in u.trajectory:
            print(ts)
            result.append(ts.positions)
        # Option to save pdbs of each trajectory
        if (save_pdbs):
            with mda.Writer("protein.pdb", protein.n_atoms) as W:
                for ts in u.trajectory:
                    W.write(protein)
    
    # Case if the file is a single pdb file
    elif file_type == "pdb":
        result.append(u.positions)

    # Return a numpy array of type float64
    result = np.array(result).astype(np.float64)
    if save_np == True:
        np.save(np_file, result)
    return result