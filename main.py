# importing libraries
import os
from data_utils import *
try:
    from pyrosetta import *
except ImportError:
    print("PyRosetta is not imported. Make sure that you have installed pyrosetta and that you are in the right environment")

# installing the required libraries
os.system('pip3 install -r requirements.txt')

# initializing PyRosetta
pyrosetta.init("-out:level 600")

# reading the PDB file
stay_gold = pyrosetta.pose_from_pdb("pdb_files/StayGold_Intein_v9.pdb")
gyra_intein = pyrosetta.pose_from_pdb("pdb_files/1am2.pdb")

# commencing analysis
# single_mutation_analysis(stay_gold, "stayGold.csv")

# single_mutation_analysis(gyra_intein, "intein.csv")
docked = pyrosetta.pose_from_pdb("pdb_files/dockedNrdj-1mCherry.pdb")
single_mutation_analysis(docked,"docked.csv")
