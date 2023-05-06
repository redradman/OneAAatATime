# importing libraries
from pyrosetta import *
from data_utils import *

# initializing PyRosetta
pyrosetta.init()

# reading the PDB file
wt_pose = pyrosetta.pose_from_pdb("StayGold_Intein_v9.pdb")


# mutant_pose = wt_pose.clone()
# print(mutant_pose)


# residue_number = 1  # Replace this with the residue number you want to mutate
# new_residue = 'ALA'    # 3 Letter code for the new residue

# mutant_pose = mutate_residue(mutant_pose, residue_number, new_residue)


# return 
# print("Wild-type score:", calculate_FA_score(wt_pose))
# print("Mutant score:", calculate_FA_score(mutant_pose))

# single_mutation_analysis(wt_pose)


# Secondary Structure Analysis code
from pyrosetta.rosetta.protocols.membrane import get_secstruct
secstruct = get_secstruct(wt_pose)

sec_strcuture = ""
for i, sec_elem in enumerate(secstruct, start=1):
    # print(f"Residue {i}: {sec_elem}")
    sec_strcuture += sec_elem


print(sec_strcuture)
    
    
    


