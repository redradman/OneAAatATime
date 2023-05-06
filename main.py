# importing libraries
from pyrosetta import *
from data_utils import *
from pyrosetta import create_score_function


# initializing PyRosetta
pyrosetta.init()

# reading the PDB file
wt_pose = pyrosetta.pose_from_pdb("StayGold_Intein_v9.pdb")


mutant_pose = wt_pose.clone()
print(mutant_pose)


# residue_number = 10  # Replace this with the residue number you want to mutate
# new_residue = 'VAL'    # 3 Letter code for the new residue

# mutant_pose = mutate_residue(mutant_pose, residue_number, new_residue)

single_mutation_analysis(wt_pose)
