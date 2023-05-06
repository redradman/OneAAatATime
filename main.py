# importing libraries
from pyrosetta import *
from data_utils import mutate_residue

# initializing PyRosetta
pyrosetta.init()

# reading the PDB file
wt_pose = pyrosetta.pose_from_pdb("StayGold_Intein_v9.pdb")


mutant_pose = wt_pose.clone()
print(mutant_pose)


residue_number = 1  # Replace this with the residue number you want to mutate
new_residue = 'ALA'    # 3 Letter code for the new residue

mutant_pose = mutate_residue(mutant_pose, residue_number, new_residue)

# Create a scoring function
score_function = get_fa_scorefxn()

# Calculate the scores
wt_score = score_function(wt_pose)
print(mutant_pose) 
mutant_score = score_function(mutant_pose)

print("Wild-type score:", wt_score)
print("Mutant score:", mutant_score)
    
    
    


