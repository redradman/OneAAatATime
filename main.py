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
    
# from pyrosetta.rosetta.core.scoring.hbonds import HBondSet

# hbond_set_wt = HBondSet(wt_pose)
# hbond_set_mutant = HBondSet(mutant_pose)

# print("Wild-type hydrogen bonds:", hbond_set_wt.nhbonds())
# print("Mutant hydrogen bonds:", hbond_set_mutant.nhbonds())



# from pyrosetta import init, pose_from_pdb
# from pyrosetta.rosetta.core.scoring import ScoreFunction, get_score_function
# from pyrosetta.rosetta.core.scoring.hbonds import HBondSet



# score_function = get_score_function()


# hbond_set = HBondSet()
# wt_pose.update_residue_neighbors()
# hbond_set.setup_for_residue_pair_energies(wt_pose, False, False)


# total_hbonds = hbond_set.nhbonds()
# print("Total hydrogen bonds:", total_hbonds)


# for i in range(1, total_hbonds + 1):
#     hbond = hbond_set.hbond(i)
#     donor_res_num = hbond.don_res()
#     acceptor_res_num = hbond.acc_res()
#     donor_atom = wt_pose.residue(donor_res_num).atom_name(hbond.don_hatm())
#     acceptor_atom = wt_pose.residue(acceptor_res_num).atom_name(hbond.acc_atm())
#     energy = hbond.energy()
#     print(f"Hydrogen bond {i}: {donor_res_num}-{donor_atom} -> {acceptor_res_num}-{acceptor_atom}, energy: {energy}")


single_mutation_analysis(wt_pose)
