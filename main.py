# importing libraries
from pyrosetta import *
from data_utils import *

# initializing PyRosetta
pyrosetta.init("-out:level 500")

# reading the PDB file
wt_pose = pyrosetta.pose_from_pdb("StayGold_Intein_v9.pdb")

# commencing analysis
single_mutation_analysis(wt_pose)
