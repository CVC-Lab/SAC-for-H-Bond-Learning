import os
from missing_residues import missing_residue
from mol_processors.pdb_processor import download_pdb

small_molecules = ['1FME', '2HBA', '2WxC', '1ENH', '1MIO', '2A3D', '1LMB', '1CQ0', '1L2Y', '1ROO', '1T5Q', '1V1D', '1LFC', '8TFV', '1KGM', '1M4F', '1TV0', '2jof', '2f4k', '1ery', '2p6j', '3gb1', '2f21', '1prb', '1bdd', '1gh1', '1fex', '1dv0']
#small_molecules = ['1PRB', '2JOF']
#TODO:Make this nicer looking
software_dir = "/home/conradli/software"
data_dir = "/home/conradli/SAC-for-H-Bond-Learning/data"
# Software paths
prepare_protein_path =  os.path.join(software_dir, "protein_prep/prepare.py")
path_to_dihe = os.path.join(software_dir,"create_psf/create_psf_dihe")
path_to_namd_conversion = os.path.join(software_dir, "charmm2namd.py")
path_to_rtf = "/home/conradli/SAC-for-H-Bond-Learning/Octree/FromBU/oct-example/params/pdbamino_new.rtf"
missing_residue_molecules = []


# Check for duplicates
for index, molecule in enumerate(small_molecules):
    for index2, molecule2 in enumerate(small_molecules):
        if molecule2.lower() == molecule.lower() and index != index2:
            print("There is duplicate of ", molecule, "at index", index)
            raise Exception("No duplicates allowed")

for molecule in small_molecules:
    molecule = molecule.upper()
    molecule_path = os.path.join(data_dir, molecule)
    if not os.path.exists(molecule_path):
        os.mkdir(molecule_path)
    
    #Input: molecule ID
    #Output: pdb in molecule's directory
    download_pdb(molecule, output_dir = molecule_path)

    #Checks for missing residues, skips iteration if there is for now
    path_to_pdb = os.path.join(molecule_path, "{}.pdb".format(molecule.lower()))
    if missing_residue(path_to_pdb):
        missing_residue_molecules.append(molecule)
        continue

    #Input: Molecule.pdb
    #Output: Molecule_pnon.pdb
    os.system("{} {}".format(prepare_protein_path, path_to_pdb))
    pnon_path =  os.path.join(molecule_path,"{}_pnon.pdb".format(molecule.lower()))
    #Input: Molecule_pnon.pdb
    #Output: Molecule_PSF
    psf_output = os.path.join(molecule_path, molecule.lower() + "_pnon.psf")
    os.system("{} {} {} {}".format(path_to_dihe, pnon_path, path_to_rtf, psf_output))

    #Input Molecule_pnon.pdb
    #Output: mol2 file
    mol2_output = os.path.join(molecule_path, molecule.lower() + "_pnon.mol2")
    os.system("obabel {} -O {}".format(pnon_path, mol2_output))
    
    # Generates NAMD verison of psf file from charmm psf file
    xplor_output = os.path.join(molecule_path, molecule.lower() + "_pnon_xplor.psf")
    os.system("python {} {} {} > {}".format(path_to_namd_conversion, psf_output, path_to_rtf, xplor_output))

print(missing_residue_molecules)
