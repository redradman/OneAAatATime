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


residue_number = 10  # Replace this with the residue number you want to mutate
new_residue = 'VAL'    # 3 Letter code for the new residue

mutant_pose = mutate_residue(mutant_pose, residue_number, new_residue)




# secstruct = get_secstruct(wt_pose)


# residue_index = 10
# print(f"Secondary structure for residue {residue_index}: {secstruct[residue_index - 1]}")

# secondary_structure = ""
# for i, sec_elem in enumerate(secstruct, start=1):
#     secondary_structure += sec_elem
    
# print("\n\n\n\n HHHHHH")
# print(secondary_structure)


print(calculate_secondary_stucture(wt_pose))

