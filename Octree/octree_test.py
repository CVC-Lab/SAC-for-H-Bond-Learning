import octree
from octree import Protein
from pympler import asizeof, tracker, refbrowser
from memory_profiler import profile
import os
# TODO
# 1. Preprocess all the pdb files
# 2. Verify that all data files work with the octree

if __name__ == '__main__':
    #@profile
    def main():
        #tr = tracker.SummaryTracker()
        print("Initialization")
        #path_prefix = "/home/conradli/SAC-for-H-Bond-Learning/Octree/FromBU/oct-example/"
        dir_prefix = "/home/jeffreymo572/Research/SAC-for-H-Bond-Learning/data"
        file_prefix = "1bdd_pnon"
        path_prefix = os.path.join(dir_prefix, "1BDD", file_prefix)
        param_prefix = os.path.join(dir_prefix, "params/")
        print(path_prefix)
        # File Paths for param and pdb files
        pdbFile = bytes(path_prefix + ".pdb", encoding="ascii")
        #pdbFixedFile = bytes(path_prefix + "-fixed.pdb", encoding="ascii") # Currently not used
        mol2File = bytes(path_prefix +".mol2", encoding="ascii")
        psfFile = bytes(path_prefix + ".psf", encoding="ascii")
        outnFile = bytes(path_prefix + "-outn.pdb", encoding="ascii") # Output for nonbonded list
        outoFile = bytes(path_prefix + "-outo.pdb", encoding="ascii") # Output for octree
        prmFile = bytes(param_prefix + "parm_new.prm", encoding="ascii")
        rtfFile = bytes(param_prefix + "pdbamino_new.rtf", encoding="ascii")
        aprmFile = bytes(param_prefix + "atoms.0.0.6.prm.ms.3cap+0.5ace.Hr0rec", encoding="ascii")
        prot = Protein (pdbFile, psfFile, mol2File, prmFile, rtfFile, aprmFile)
        print(prot.extract_features())
        print("After min")
    #print(octree.octree_test())
    main()

