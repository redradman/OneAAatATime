from pyrosetta import *
# from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
# from pyrosetta.toolbox.mutants import mutate_residue
from pyrosetta.toolbox import mutate_residue

pyrosetta.init()


wt_pose = pyrosetta.pose_from_pdb("StayGold_Intein_v9.pdb")
# print(wt_pose)


mutant_pose = wt_pose.clone()

residue_number = 11  # Replace this with the residue number you want to mutate
new_residue = 'B'    # Replace 'A' with the one-letter code of the desired amino acid


mutant_pose = mutate_residue(mutant_pose, residue_number, new_residue)

# Create a scoring function
score_function = get_fa_scorefxn()

# Calculate the scores
wt_score = score_function(wt_pose)
print(mutant_pose) # ERROR
mutant_score = score_function(mutant_pose)

print("Wild-type score:", wt_score)
print("Mutant score:", mutant_score)
