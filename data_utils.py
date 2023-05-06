# importing libraries 
import pandas as pd
from pyrosetta import *
from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
from pyrosetta.rosetta.protocols.membrane import get_secstruct
from pyrosetta.rosetta.core.scoring.hbonds import HBondSet
from pyrosetta.rosetta.core.scoring import calc_total_sasa


# added for gyration 
# import math
# from pyrosetta.rosetta.numeric import xyzVector_double_t



# initializing PyRosetta
pyrosetta.init()


################################################################
def single_mutation_analysis(wt_pose):
    """
    Creates all possible single AA mutation when given a wildtype pose
    """
    
    
    # Creating the pandas data structure
    df = pd.DataFrame(columns=['type', # wild_type or mutant
                               'residue_number', # the number of the targeted residue
                               'previous_aa', # Previous amino acid (in the wildtype)
                               'new_aa_1l', # the 1 letter representation of the new amino acid
                               'new_aa_3l', # the 3 letter representation of the new amino acid
                               'conversion', # resiude number + previous aa + new aa i.e. 3AtoR (at resiude number 3 A was converted to R)
                               'new_seq', # the new amino acid sequence after the mutation has been applied (if it is mutant otherwise it is wildtype)
                               'fa_score', # the Full Atom score for the new amino acid
                               'ddg_score', # the ddG (change in binding free energy) score
                               'hbond_score', # hydrogen bonding patterns
                               'sasa_score', # solvent accessible surface area
                               'secondary_structure'
                               ])
    
    # adding the wildtype in the 1st row of the data frame
    df.loc[len(df)] = ['wild_type', 'NA', 'NA', 'NA',  'NA',  'NA', wt_pose.sequence(), calculate_FA_score(wt_pose),
                               0, # the ddG (change in binding free energy)
                               calculate_hbonds_comprehensive(wt_pose), # hydrogen bonding patterns
                               calc_sasa(wt_pose), # solvent accessible surface area
                               'secondary_structure'
                               ]
    
    
    
    # converting the data frame to a csv file
    df.to_csv('example.csv', index=False)


################################################################
def calculate_FA_score(pose):
    """ 
    Returns the Full Atom score for the passed pose.
    """
    # Create a scoring function
    score_function = get_fa_scorefxn() 
    # score_function = create_score_function("score12")
    return score_function(pose)


################################################################
def calculate_ddg_score(wt_pose, mut_pose):
    """
    Returns the ddG (change in binding free energy) score for the passed pose.
    """
    # Create a scoring function
    score_function = create_score_function('ref2015') 
    
    # Calculate score the wild-type and mutant poses
    wt_score = score_function(wt_pose)
    mut_score = score_function(mut_pose)
    
    # calculate ddg
    ddg_score = mut_score - wt_score
    
    return ddg_score

################################################################
def calculate_secondary_stucture(pose) -> str:
    """
    returns the secondary structure of the passed pose in string format
    """
    secstruct = get_secstruct(pose)
    sec_strcuture = ""
    for i, sec_elem in enumerate(secstruct, start=1):
        sec_strcuture += sec_elem
    return sec_strcuture

################################################################
def calculate_hbonds_simple(pose) -> int:
    """
    Returns the number of hbonds in the passed pose
    this function calculates the number of hbonds in a simple manner and is computationally less expensive
    
    calculates hydrogen bonds based on distance and angle criteria
    This method is suitable for analyzing small to medium-sized proteins 
    or for screening large numbers of protein structures for hydrogen bonding.
    """
    hbond_set_wt = HBondSet(pose)
    return hbond_set_wt.nhbonds()


################################################################
def calculate_hbonds_comprehensive(pose) -> int:
    """
    Returns the number of hbonds in the passed pose
    this function calculates the number of hbonds in a more detailed manner
    
    This method takes into account the orientation and position of each residue in the protein,
    as well as the potential for multiple hydrogen bonds to form between different residues.
    This method is suitable for analyzing large or complex protein structures,
    or for performing detailed analyses of specific hydrogen bonding interactions
    """
    hbond_set = HBondSet()
    pose.update_residue_neighbors()
    hbond_set.setup_for_residue_pair_energies(pose, False, False)


    total_hbonds = hbond_set.nhbonds()
    return total_hbonds


################################################################
def iterate_through_all_hbonds(pose):
    """
    prints all of the hbonfs in the passed pose
    """
    score_function = get_score_function()

    hbond_set = HBondSet()
    pose.update_residue_neighbors()
    hbond_set.setup_for_residue_pair_energies(pose, False, False)

    total_hbonds = hbond_set.nhbonds()
    print("Total hydrogen bonds:", total_hbonds)

    for i in range(1, total_hbonds + 1):
        hbond = hbond_set.hbond(i)
        donor_res_num = hbond.don_res()
        acceptor_res_num = hbond.acc_res()
        donor_atom = pose.residue(donor_res_num).atom_name(hbond.don_hatm())
        acceptor_atom = pose.residue(acceptor_res_num).atom_name(hbond.acc_atm())
        energy = hbond.energy()
        print(f"Hydrogen bond {i}: {donor_res_num}-{donor_atom} -> {acceptor_res_num}-{acceptor_atom}, energy: {energy}")

################################################################
######### FUNCTION BRLOW IS NOT WORKING 
# def calculate_gyration_radius(pose) -> float:
#     """
#     Returns the radius of gyration for the passed pose.
    
#     NOT FUNCTIONING
#     """
    
#     def center_of_mass(pose):
#         total_mass = 0.0
#         com = xyzVector_double_t(0.0, 0.0, 0.0)

#         for i in range(1, pose.total_residue() + 1):
#             residue = pose.residue(i)
#             for j in range(1, residue.natoms() + 1):
#                 atom = residue.atom(j)
#                 mass = residue.atom_type_set().atom_type(residue.atom_type_index(j)).mass()
#                 com += atom.xyz() * mass
#                 total_mass += mass

#         return com / total_mass
    
#     def radius_of_gyration(pose):
#         com = center_of_mass(pose)
#         num_residues = pose.total_residue()
#         rg = 0.0

#         for i in range(1, num_residues + 1):
#             residue = pose.residue(i)
#             for j in range(1, residue.natoms() + 1):
#                 atom = residue.atom(j)
#                 rg += atom.xyz().distance_squared(com)

#         rg = math.sqrt(rg / pose.total_atoms())
#         return rg
    
#     gyration_radius = radius_of_gyration(wt_pose)

#     return gyration_radius



################################################################    
def calc_sasa(pose) -> float:
    """ 
    calculates the solvent accessible surface area of a pose for water
    """
    probe_radius = 1.4  # Use the default probe_radius of 1.4 Å for water
    sasa = calc_total_sasa(pose, probe_radius)
    return sasa
    
    
def calc_sasa(pose, probe_radius) ->float:
    """
    calculates the solvent accessible surface area for a pose given a radius 
    The probe radius represents the size of the solvent molecule
    
    The probe_radius is a parameter used in the calculation of the solvent accessible surface area (SASA) of a biomolecule,
    such as a protein. 
    The SASA calculation method is based on the rolling sphere algorithm, 
    where a spherical "probe" of a certain radius rolls over the surface of the molecule.
    """
    return calc_total_sasa(pose, probe_radius)

    


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
