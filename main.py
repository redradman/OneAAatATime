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


# return 
# print("Wild-type score:", calculate_FA_score(wt_pose))
# print("Mutant score:", calculate_FA_score(mutant_pose))


# single_mutation_analysis(wt_pose)
# iterate_through_all_hbonds(wt_pose)

# from pyrosetta import create_score_function

# # Create a score function
# score_function = create_score_function("ref2015")


# # Score the wild-type and mutant poses
# wt_score = score_function(wt_pose)
# mut_score = score_function(mutant_pose)

# # Calculate the ∆∆G score
# ddg_score = mut_score - wt_score

# # Print the calculated ∆∆G score
# print("∆∆G score:", ddg_score)

# # print(calculate_FA_score(wt_pose))
# # print(calculate_FA_score(mutant_pose))
# print(calculate_FA_score(mutant_pose) - calculate_FA_score(wt_pose))

# score_function1 = get_fa_scorefxn() 
# score_function2 = create_score_function("score12")
# print(score_function1(wt_pose))
# print(score_function2(wt_pose))

from pyrosetta.rosetta.core.scoring import calc_total_sasa


# Define SASA options


# Calculate SASA
wt_sasa = calc_total_sasa(wt_pose, 1.4)
mutant_sasa = calc_total_sasa(mutant_pose, 1.4)

print("Wild-type SASA:", wt_sasa)
print("Mutant SASA:", mutant_sasa)




