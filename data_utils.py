# importing libraries 
from pyrosetta import *
from pyrosetta.rosetta.protocols.simple_moves import MutateResidue

# initializing PyRosetta
pyrosetta.init()


################################################################









################################################################
def mutate_residue(pose, residue_number, new_residue):
    """ 
    creates a new pose with a the given residue number and a new residue.
    
    For example, mutate_residue("AAA", 1, R) return the new pose "RAA"
    """
    mutate = MutateResidue(target=residue_number, new_res=new_residue)
    mutate.apply(pose)
    return pose



################################################################
def find_unique_aa_in_pdb(pose):
    """
    Returns a set of unique amino acids in the passed PDB file
    """
    unique_amino_acids = set()
    
    for i in range(1, pose.total_residue() + 1):
        res = pose.residue(i)
        if res.is_protein():
            unique_amino_acids.add(res.name3())
            
    return f"Unique amino acids in the PDB file: {unique_amino_acids}"


# print(mutant_pose.sequence())
