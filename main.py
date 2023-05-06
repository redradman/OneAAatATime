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


# import math
# from pyrosetta.rosetta.numeric import xyzVector_double_t



# def center_of_mass(pose):
#     total_mass = 0.0
#     com = xyzVector_double_t(0.0, 0.0, 0.0)

#     for i in range(1, pose.total_residue() + 1):
#         residue = pose.residue(i)
#         for j in range(1, residue.natoms() + 1):
#             atom = residue.atom(j)
#             mass = residue.atom_type_set().atom_type(residue.atom_type_index(j)).mass()
#             com += atom.xyz() * mass
#             total_mass += mass

#     return com / total_mass

# def radius_of_gyration(pose):
#     com = center_of_mass(pose)
#     num_residues = pose.total_residue()
#     rg = 0.0

#     for i in range(1, num_residues + 1):
#         residue = pose.residue(i)
#         for j in range(1, residue.natoms() + 1):
#             atom = residue.atom(j)
#             rg += atom.xyz().distance_squared(com)

#     rg = math.sqrt(rg / pose.total_atoms())
#     return rg



# gyration_radius = radius_of_gyration(wt_pose)

# print("Radius of gyration:", gyration_radius)



def get_atom_mass(pose, residue_index, atom_index):
    residue = pose.residue(residue_index)
    atom_type_index = residue.atom_type_index(atom_index)
    atom_type_set = pose.conformation().atom_type_set()
    atom_type = atom_type_set[atom_type_index]
    mass = atom_type.mass()
    return mass


residue_index = 1
atom_index = 1
mass = get_atom_mass(wt_pose, residue_index, atom_index)

print("Mass of atom at residue index", residue_index, "and atom index", atom_index, "is:", mass)

