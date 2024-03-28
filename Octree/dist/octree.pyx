import numpy as np
cimport numpy as np
from libc.stdlib cimport malloc, free
from libc.stdint cimport uintptr_t


np.import_array()

cdef extern from "numpy/arrayobject.h":
    void PyArray_ENABLEFLAGS(np.ndarray arr, int flags)
    
# Import functions and structs from octree c code
cdef extern from "mol.0.0.6.h":
    cdef struct atombond:
        int ai
        int aj
    ctypedef atombond mol_bond
    cdef struct atomangle:
        float th
    cdef struct atomtorsion:
        float chi
    cdef struct atomimproper:
        float psi
    cdef struct atom:
        int num_neighbors
        int *neighbors
        int nbondis
        int *bondis
        int nangs
        char *name
        int ingrp
        atomangle ** angs
        int ntors
        atomtorsion ** tors
        int nimps
        atomimproper ** imps
        double chrg
        double rminh
        double X
        double Y
        double Z
    cdef struct atomgrp:
        int nactives
        int natoms
        atom * atoms
        int nbonds
        mol_bond bonds
    cdef struct prmatom:
        char *typemaj
        char *typemin
        float r
        float q
    cdef struct prm:
        int natoms
        prmatom * atoms
    cdef struct my_par
    cdef void ag2array( double* array, atomgrp * ag)
    cdef void array2ag ( double* array, atomgrp * ag)
    cdef atomgrp * create_atomgrp (char *pdbFile, char *pdbFixedFile, char *psfFile, char *mol2File, char *prmFile, char *rtfFile, char *aprmFile, int useHbond) 
    cdef void hydro_min (atomgrp * ag, char* pdbFile, char *outnFile, char *outoFile, int maxIter, int useHbond, double *dcHbond, double *dcVdw, double *dcElec)
    cdef void my_en_grad( int n, double *inp, void *prms, double *en, double *grad )
    cdef my_par * create_par(atomgrp * ag, int useVdw, int useElec, int useHbond ,double *dcVdw, double *dcElec, double *dcHbond, double aprxVdw, double aprxElec)
    cdef int * get_bondls (atomgrp * ag, atom * a)
    cdef void only_find_neighbors(int n, double *inp, void *prms)

cdef class Protein:
    cdef atomgrp * ag_ptr
    cdef my_par * my_parms
    cdef double * atom_coords
    cdef dict __dict__
    cdef char *pdbFile
    cdef char *psfFile 
    cdef char *mol2File
    cdef char *prmFile 
    cdef char *rtfFile 
    cdef char *aprmFile
    cdef int neighbor_type
    cdef double dcVdw[1]
    cdef double dcElec[1]
    cdef double dcHbond[1]
    cdef double energy[1]
    cdef int maxIter
    cdef int noOct
    cdef int noNblist
    cdef int useVdw
    cdef int useElec
    cdef int useHbond
    cdef double aprxVdw
    cdef double aprxElec

    def __init__(self, char *pdbFile, char *psfFile, char *mol2File, char *prmFile, char *rtfFile, char *aprmFile):
        # Set Params
        path_prefix = "/home/jeffreymo572/Research/SAC-for-H-Bond-Learning/Octree/FromBU/oct-example/"
        self.maxIter = 100
        #self.noOct = 0
        #self.noNblist = 0
        self.useVdw = 0
        self.useElec = 0
        self.useHbond = 1
        self.aprxVdw = 1.0
        self.aprxElec = 1.0
        self.dcVdw[0] = 12.000000
        self.dcElec[0] = 15.000000
        self.dcHbond[0] = 3.00000
        self.energy[0] = 0 
        self.neighbor_type = 1
        self.pdbFile = pdbFile
        self.psfFile = psfFile
        self.mol2File = mol2File
        self.prmFile = prmFile
        self.rtfFile = rtfFile
        self.aprmFile = aprmFile
        #self.outnFile = outnFile
        #self.outoFile = outoFile

        # Create an atom group from the given inputs
        self.ag_ptr = create_atomgrp (self.pdbFile, NULL, self.psfFile, self.mol2File, self.prmFile, self.rtfFile, self.aprmFile, self.useHbond)
        # Create octree parms and bulild octree
        self.my_parms = create_par(self.ag_ptr, self.useVdw, self.useElec, self.useHbond, self.dcVdw, self.dcElec, self.dcHbond, self.aprxVdw, self.aprxElec)
       
        # Convert atom group coordinates to an array of doubles
        self.num_atoms = self.ag_ptr.natoms
        self.active_atoms_dims = 3 * self.ag_ptr.natoms
        self.atom_coords = <double *> malloc(self.active_atoms_dims * sizeof(double)) 
        
        ag2array(self.atom_coords, self.ag_ptr)
        # Calculate the initial energy and fill neighbor lists
        my_en_grad(self.active_atoms_dims, self.atom_coords, self.my_parms, self.energy, NULL)
        # Grab coords as numpy array
        self.cur_coords = ptr_to_nparray_double(self.atom_coords, self.active_atoms_dims)
        #self.update_energy (self.cur_coords)
        # Grab neighborhood lists
        self.neighbor_lists, self.num_neighbor_lists = get_neighbor_lists (self.ag_ptr, self.neighbor_type)
        
        # USED FOR DEBUGGING
        #hydro_min (self.ag_ptr, self.pdbFile, self.outnFile, self.outoFile, self.maxIter, self.useHbond, self.dcHbond, self.dcVdw, self.dcElec)
        #self.validate_neighbors()
        #print("This is the energy", self.energy)
    
    # Only updates the neighbor lists
    def only_update_neighbors(self, new_coords, neigh_type=1):
        only_get_neighbors(new_coords, self.active_atoms_dims, self.my_parms)
        self.neighbor_lists, self.num_neighbor_lists = get_neighbor_lists (self.ag_ptr, neigh_type)
        return self.neighbor_lists, self.num_neighbor_lists

    # Recalculates the energy based on the new coordinates
    # Params:
    #
    # new_coords: a np.float64 numpy array of size (num_atoms, 3) 
    # Returns the energy, a python array of neighbor lists, and a python array containing the number of neighbors in each neighbor list
    def update_energy(self, new_coords):
        calc_energy (new_coords, self.active_atoms_dims, self.my_parms, self.energy)
        self.neighbor_lists, self.num_neighbor_lists = get_neighbor_lists (self.ag_ptr, self.neighbor_type)
        return self.energy[0], self.neighbor_lists, self.num_neighbor_lists
    
    # Updates cur_coords. We do this because we have to 
    def update_cur_coords(self):
        self.cur_coords = ptr_to_nparray_double(self.atom_coords, self.active_atoms_dims)
    
    def extract_features(self):
        return get_atomistic_features(self.ag_ptr)

# Only finds neighbors
cdef only_get_neighbors(np.float64_t[:] new_coords, int active_atoms_dims, void *prms):
    cdef np.float64_t* addr = &new_coords[0]
    only_find_neighbors(active_atoms_dims, addr, prms)
    return

# Calculates the energy using the new coordinates from a numpy array of coordinates
cdef calc_energy (np.float64_t[:] new_coords, int active_atoms_dims, void *prms, double *en):
    cdef np.float64_t* addr = &new_coords[0]
    my_en_grad(active_atoms_dims, addr, prms, en, NULL)
    return

cdef get_atomistic_features(atomgrp * ag_ptr):
    charges = [None] * ag_ptr.natoms
    radii = [None] * ag_ptr.natoms
    names = [None] * ag_ptr.natoms
    for i in range(ag_ptr.natoms):
        charge = (&ag_ptr.atoms[i]).chrg
        vdw_radius = (&ag_ptr.atoms[i]).rminh
        atom_name = (&ag_ptr.atoms[i]).name[0]
        names[i] = atom_name
        charges[i] = charge
        radii[i] = vdw_radius
    return names, charges, radii

# Gathers all atom neighborhoods into one list
# Neighborhoods can either be defined by the h-bond dist cutoff or the bond structure
# Params:
# ag_ptr: atomgroup pointer
# neighbor_type: type of neighborhood
#   1: base neighbors on dist cutoff
#   anything else: base neighbors on covalent bonds 
cdef get_neighbor_lists (atomgrp * ag_ptr, int neighbor_type):
    neighbor_lists = [None] * ag_ptr.natoms
    num_neighbors_lists = [None] * ag_ptr.natoms
    # Our neighborhood is based on the hbond dist cutoff
    if neighbor_type == 1:
        for i in range (ag_ptr.natoms):
            num_neighbors = (&ag_ptr.atoms[i]).num_neighbors
            neighbors = get_distcut_neighbors(&ag_ptr.atoms[i])
            neighbor_lists[i] = neighbors
            num_neighbors_lists[i] = num_neighbors
    # Our neighborhood is based on the bond structure
    else:
        for i in range (ag_ptr.natoms):
            num_bonds = (&ag_ptr.atoms[i]).nbondis
            bonds = ptr_to_nparray_int(get_bondls(ag_ptr, &ag_ptr.atoms[i]), num_bonds)
            neighbor_lists[i] = bonds
            num_neighbors_lists[i] = num_bonds
    return neighbor_lists, num_neighbors_lists
    
# Function to check that neighbors are or are not in each other's lists
def check_nblist_doubled(neighbor_lists, doubled=False):
    for atom_index, array in enumerate(neighbor_lists):
        for neighbor_index in array:
            # Check if the atom which owns this neighborhood list is in its neighbor's neighborhood list
            if doubled:
                if (atom_index not in neighbor_lists[neighbor_index]):
                    print (neighbor_index, " and ", atom_index, "do not contain each other")
                    return False
            # Otherwise check the oppposite
            else:
                if (atom_index in neighbor_lists[neighbor_index]):
                    print (neighbor_index, " and ", atom_index, "contain each other")
                    return False
    return True

# Grabs the neighborhood list of atom
# Returns an empty list if the atom has no neighbors
#TODO there could be a memory leak here
cdef get_distcut_neighbors(atom * a):
    
    if (a.num_neighbors == 0):
        return np.array([])
    else:
        np_a = np.copy(ptr_to_nparray_int (a.neighbors, a.num_neighbors))
        #print(np_a.flags)
        #np_a.free_data = True
        return np_a

''' ############################Helper Utilities############################ '''


# Checks for duplicates in an array
# Returns if the array has duplicates
def has_duplicates(array):
    duplicate = len(array) != len(set(array))
    if (duplicate):
        print ("There is a duplicate")
    return duplicate

# Converts a double pointer to an array
# O(1) time complexity and in place
cdef ptr_to_nparray_double(double * ptr, int size):
    if size == 0:
        raise RuntimeError('size cannot be 0')
    cdef double[:] view = <double[:size]> ptr
    return np.asarray(view)

# Converst int pointer to an array
# O(1) time complexity and in place
cdef ptr_to_nparray_int(int * ptr, int size):
    if size == 0:
        raise RuntimeError('size cannot be 0')
    cdef int[:] view = <int[:size]> ptr
    return np.asarray(view)

# Converts double array to python list
# O(N) time complexity and copies to new array
cdef doublearray_to_python (double *ptr, int length):
    cdef int i
    lst = []
    # Copy the results into a python list
    for i in range(length):
        lst.append (ptr[i])
    # Free the memory
    if (ptr is not NULL):
        free (ptr)
    return lst

# Converts integer array to python list
# O(N) time complexity and cipies to new array
cdef intarray_to_python (int *ptr, int length):
    cdef int i
    lst = []
    # Copy the results into a python list
    for i in range(length):
        lst.append (ptr[i])
    # Free the memory
    if (ptr is not NULL):
        free (ptr)
    return lst
