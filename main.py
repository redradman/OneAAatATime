# importing libraries
from pyrosetta import *
from data_utils import *

# initializing PyRosetta
pyrosetta.init("-out:level 600")

# reading the PDB file
stay_gold = pyrosetta.pose_from_pdb("StayGold_Intein_v9.pdb")
gyra_intein = pyrosetta.pose_from_pdb("1am2.pdb")

# commencing analysis
# single_mutation_analysis(stay_gold, "stayGold.csv")

single_mutation_analysis(gyra_intein, "gyra_intein.csv")
