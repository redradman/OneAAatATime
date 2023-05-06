from pyrosetta import *
from pyrosetta.rosetta.protocols.simple_moves import MutateResidue

pyrosetta.init()


wt_pose = pyrosetta.pose_from_pdb("StayGold_Intein_v9.pdb")


mutant_pose = wt_pose.clone()
print(mutant_pose)


residue_number = 1  # Replace this with the residue number you want to mutate
new_residue = 'ALA'    # 3 Letter code for the new residue


def mutate_residue(pose, residue_number, new_residue):
    mutate = MutateResidue(target=residue_number, new_res=new_residue)
    mutate.apply(pose)
    return pose

mutant_pose = mutate_residue(mutant_pose, residue_number, new_residue)

# Create a scoring function
score_function = get_fa_scorefxn()

# Calculate the scores
wt_score = score_function(wt_pose)
print(mutant_pose) # ERROR
mutant_score = score_function(mutant_pose)

print("Wild-type score:", wt_score)
print("Mutant score:", mutant_score)
