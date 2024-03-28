import os 

run_charmm_origin_path = "/home/jeffreymo572/Research/SAC-for-H-Bond-Learning/data/1BDD/run_charmm.inp"
path_to_charmm_app = "/home/jeffreymo572/charmm/charmm/bin/charmm"
molecule_list = ["1BDD"]


for molecule in molecule_list:
    molecule = molecule.lower()
    charmm_file = open(run_charmm_origin_path, "r")
    list_of_lines = charmm_file.readlines()

    #Line 9 is .pdb input name
    list_of_lines[9] = f"open unit 1 card read name {molecule}.pdb\n"

    #Line 36 is .pdb output name
    list_of_lines[36] = f"write coor pdb name {molecule}_pnon_charmm.pdb\n"

    #Line 40 is .psf output name
    list_of_lines[40] = f"write psf CARD OLDPSF name {molecule}_pnon_charmm.psf\n"

    #Line 44 is .crd output name
    list_of_lines[44] = f"write coor card name {molecule}_pnon_charmm.crd"

    charmm_file = open(run_charmm_origin_path, "w")
    charmm_file.writelines(list_of_lines)
    charmm_file.close()

    os.system(f"{path_to_charmm_app} < {run_charmm_origin_path}")


