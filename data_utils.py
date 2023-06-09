# importing libraries 
import pandas as pd
from pyrosetta import *
from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
from pyrosetta.rosetta.protocols.membrane import get_secstruct
from pyrosetta.rosetta.core.scoring.hbonds import HBondSet
from pyrosetta.rosetta.core.scoring import calc_total_sasa
from pyrosetta.rosetta.protocols.membrane import get_secstruct
from tqdm import tqdm


from pyrosetta.rosetta.core.conformation import ResidueFactory
from pyrosetta.rosetta.core.chemical import ChemicalManager


from pyrosetta import rosetta as core

amino_acids = "ACDEFGHIKLMNPQRSTVWY"  # string of 20 standard amino acids

# added for gyration 
# import math
# from pyrosetta.rosetta.numeric import xyzVector_double_t
################################################################
def init(wt_pose):
    """Creates the data frame and the calculates the scores for wild_type"""
    # global df 
    # df = make_data_frame()
    global wt_hbonds
    global wt_sasa 
    global wt_secondary
    wt_hbonds = calculate_hbonds_simple(wt_pose)
    wt_sasa = calc_sasa_water(wt_pose)
    wt_secondary = calculate_secondary_stucture(wt_pose)

################################################################
def single_mutation_analysis(wt_pose, filename: str):
    """
    Creates all possible single AA mutation when given a wildtype pose
    """
    # setting up
    init(wt_pose)
    
    # make df and add wild type
    df = make_data_frame()
    add_wildtype(wt_pose, df)
    # the loop for adding all of the mutants

    for i in tqdm(range(1, wt_pose.total_residue() + 1), desc="Mutating residues"):
    # for i in range(1, wt_pose.total_residue() + 1):
        current_residue = wt_pose.residue(i).name1()
        for aa in amino_acids:
            if aa != current_residue:
                mutant_pose = wt_pose.clone()
                residue = wt_pose.residue(i)
                previous_aa = residue.name1()
                mutant_pose = mutate_residue(mutant_pose, i, get_three_letter_code(aa))
                
                mut_hbonds = calculate_hbonds_simple(mutant_pose)
                mut_sasa = calc_sasa_water(mutant_pose)
                mut_secondary = calculate_secondary_stucture(mutant_pose)
                
                # adding row to the data frame
                df.loc[len(df)] = ['point_mutantion', # wild_type or mutant type
                               i,  # residue
                               previous_aa, # previous aa
                               aa, # new AA 1L
                               get_three_letter_code(aa), # new AA 3L
                               str(i) + previous_aa + "to" + aa, # resiude number + previous aa + new aa i.e. 3AtoR (at resiude number 3 A was converted to R)
                               mutant_pose.sequence(), # new seq
                               calculate_FA_score(mutant_pose), # FA score
                               calculate_ddg_score(wt_pose, mutant_pose), # the ddG score
                               mut_hbonds, # hydrogen bonding
                               mut_sasa, # SASA
                               mut_secondary, # 2nd-ary structure string
                               mut_hbonds - wt_hbonds, # diff hbonds
                               mut_sasa - wt_sasa, # diff sasa
                               string_difference(wt_secondary, mut_secondary) # diff 2nd-ary structure
                               ]
                
    
    # converting the data frame to a csv file
    df.to_csv(filename, index=False)
    
def single_insertion(wt_pose, filename: str):
    """ 
    Calculates the score for all of the possible single insertion mutations of the wild-type amino acid sequence of the pdb file
    """
    # setting up
    init(wt_pose)
    
    # make df and add wild_type
    df = make_data_frame()
    add_wildtype(wt_pose, df)

    # get wild type sequence and make a list of all of the possible mutants sequences
    wt_seq = wt_pose.sequence()
    single_insertion_seqs = single_insertion_sequences(wt_seq)
    
    for i in tqdm(single_insertion_seqs, desc="Mutating residues"):
        residue =i["insert_position"]
        previous_aa = i["previous_aa"]
        new_aa = i['new_aa']
        mutant_pose = pose_from_sequence(i["new_sequence"], auto_termini=True)
        
        mut_hbonds = calculate_hbonds_simple(mutant_pose)
        mut_sasa = calc_sasa_water(mutant_pose)
        mut_secondary = calculate_secondary_stucture(mutant_pose)
        
        # adding row to the data frame
        df.loc[len(df)] = ['single_aa_insertion', # wild_type or mutant type
            residue,  # residue
            "NA", # previous aa
            new_aa, # new AA 1L
            get_three_letter_code(new_aa), # new AA 3L
            "NA", # resiude number + previous aa + new aa i.e. 3AtoR (at resiude number 3 A was converted to R)
            mutant_pose.sequence(), # new seq
            calculate_FA_score(mutant_pose), # FA score
            calculate_ddg_score(wt_pose, mutant_pose), # the ddG score
            mut_hbonds, # hydrogen bonding
            mut_sasa, # SASA
            mut_secondary, # 2nd-ary structure string
            mut_hbonds - wt_hbonds, # diff hbonds
            mut_sasa - wt_sasa, # diff sasa
            string_difference(wt_secondary, mut_secondary) # diff 2nd-ary structure
            ]
        
    df.to_csv(filename, index=False)
    

def single_deletion(wt_pose, filename: str):
    """
    Creates all of the possible amino acid sequences where 1 aa is removed from wild type and calculates the score for each of them
    """
    # setting up
    init(wt_pose)
    
    # make df and add wild_type
    df = make_data_frame()
    add_wildtype(wt_pose, df)
    
    # get wild type sequence and make a list of all of the possible mutants sequences
    wt_seq = wt_pose.sequence()
    single_insertion_seqs = single_deletion_sequences(wt_seq)
    
    for i in tqdm(single_insertion_seqs, desc="Mutating residues"):
        residue =i["deletion_position"]
        previous_aa = i["previous_aa"]
        mutant_pose = pose_from_sequence(i["new_sequence"], auto_termini=True)
        
        mut_hbonds = calculate_hbonds_simple(mutant_pose)
        mut_sasa = calc_sasa_water(mutant_pose)
        mut_secondary = calculate_secondary_stucture(mutant_pose)
        
        # adding row to the data frame
        df.loc[len(df)] = ['single_aa_deletion', # wild_type or mutant type
            residue,  # residue
            previous_aa, # previous aa
            "NA", # new AA 1L
            "NA", # new AA 3L
            "NA", # resiude number + previous aa + new aa i.e. 3AtoR (at resiude number 3 A was converted to R)
            mutant_pose.sequence(), # new seq
            calculate_FA_score(mutant_pose), # FA score
            calculate_ddg_score(wt_pose, mutant_pose), # the ddG score
            mut_hbonds, # hydrogen bonding
            mut_sasa, # SASA
            mut_secondary, # 2nd-ary structure string
            mut_hbonds - wt_hbonds, # diff hbonds
            mut_sasa - wt_sasa, # diff sasa
            string_difference(wt_secondary, mut_secondary) # diff 2nd-ary structure
            ]
        
    df.to_csv(filename, index=False)
    
    

def single_insertion_sequences(sequence):
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"  # string of 20 standard amino acids
    new_sequences = []

    for i in range(len(sequence) + 1):
        for amino_acid in amino_acids:
            new_sequence = sequence[:i] + amino_acid + sequence[i:]
            previous_char = sequence[i-1] if i > 0 else sequence[-1]
            new_sequences.append({
                'new_sequence': new_sequence,
                'previous_aa': previous_char,
                'insert_position': i + 1,
                'new_aa': amino_acid
            })

    return new_sequences

def single_deletion_sequences(sequence):
    """ 
    Creates all of the mutations of the passed sequence such that in each iteration a single amino acid is deleted from the squence.
    """
    new_sequences = []

    for i in range(len(sequence)):
        new_sequence = sequence[:i] + sequence[i+1:]
        removed_char = sequence[i]
        new_sequences.append({
            'new_sequence': new_sequence,
            'previous_aa': removed_char,
            'deletion_position': i + 1
        })
    print(new_sequences)
    return new_sequences

def make_data_frame():
    """ 
    creates a pandas data frame which will be used as a skeleton to store all of the data
    """
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
                               'secondary_structure',
                               'diff_hbonds',
                               'diff_sasa',
                               'diff_secondary_structure'
                               ])
    
    return df
    
def add_wildtype(wt_pose, df):
    """ 
    adds the corresponding values of wildtype to the csv file
    """
    wt_seq = wt_pose.sequence()
    wt_FA_score = calculate_FA_score(wt_pose)
    wt_hbonds = calculate_hbonds_simple(wt_pose)
    wt_sasa = calc_sasa_water(wt_pose)
    wt_secondary = calculate_secondary_stucture(wt_pose)
    
    # adding the wildtype in the 1st row of the data frame
    df.loc[len(df)] = ['wild_type', 
                       'NA', 
                       'NA', 
                       'NA',  
                       'NA',  
                       'NA', 
                       wt_seq, 
                       wt_FA_score,
                        0, # the ddG 
                        wt_hbonds, # hydrogen bonding 
                        wt_sasa, # SASA
                        wt_secondary,
                        0,
                        0,
                        0
                        ]
    

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
def calc_sasa_water(pose) -> float:
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
def calculate_secondary_stucture(pose) -> str:
    """
    returns the secondary structure of the passed pose in string format
    """
    secstruct = get_secstruct(pose)
    sec_strcuture = ""
    for _, sec_elem in enumerate(secstruct, start=1):
        sec_strcuture += sec_elem
    return sec_strcuture



def string_difference(str1, str2) -> int:
    """
    Computes the difference between two strings by counting the number of characters that differ between them.
    Args:

    Returns:
        int: The number of characters that differ between the two strings.
    """
    # Ensure that the input strings have the same length
    # if len(str1) != len(str2):
    #     raise ValueError("Input strings must have the same length.")
    
    # Check if the strings are the same
    if str1 == str2:
        return 0
    # Compute the difference by counting the number of differing characters
    diff_count = sum(1 for c1, c2 in zip(str1, str2) if c1 != c2)
    return diff_count


################################################################
def mutate_residue(pose, residue_number, new_residue):
    """ 
    creates a new pose with a the given residue number and a new residue.
    
    For example, mutate_residue("AAA", 1, VAL) return the new pose "VAA"
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
