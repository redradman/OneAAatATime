# importing libraries 
from pyrosetta import *
from pyrosetta.rosetta.protocols.simple_moves import MutateResidue

# initializing PyRosetta
pyrosetta.init()


################################################################
def single_mutation_analysis(wt_pose):
    """
    Creates all possible single AA mutation when given a wildtype pose
    """
    pass


################################################################
def calculate_FA_score(pose):
    """ 
    Returns the Full Atom score for the passed pose.
    """
    # Create a scoring function
    score_function = get_fa_scorefxn() 
    return score_function(pose)


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



################################################################
def get_three_letter_code(one_letter_code: str) -> str:
    """
    Returns the 3 letter code for the passed amino acid 1 letter code
    """
    aa_codes = {
        'A': 'ALA',
        'C': 'CYS',
        'D': 'ASP',
        'E': 'GLU',
        'F': 'PHE',
        'G': 'GLY',
        'H': 'HIS',
        'I': 'ILE',
        'K': 'LYS',
        'L': 'LEU',
        'M': 'MET',
        'N': 'ASN',
        'P': 'PRO',
        'Q': 'GLN',
        'R': 'ARG',
        'S': 'SER',
        'T': 'THR',
        'V': 'VAL',
        'W': 'TRP',
        'Y': 'TYR'
    }

    return aa_codes.get(one_letter_code.upper(), None)

# print(mutant_pose.sequence())
