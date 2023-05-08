# OneAAatATime
Before running the `main.py` make sure that you have installed [PyRosetta](https://www.pyrosetta.org/downloads) and [python3](https://www.python.org/downloads/).
You can clone this repo and add your .pdb file to the same directory as `main.py` to simulate single AA mutation across all of the AA chain of your .pdb file. 

## How to do single amino acid mutation on a new .pdb file
**Add the lines below to `main.py` to simulate the single amino acid mutation across the entire chain.**
```python
new_protein = pyrosetta.pose_from_pdb("new_protein.pdb")

single_mutation_analysis(new_protein, "data.csv")
```
`single_mutation_analysis()` takes **pose** and **filename**. Afterwards, each amino acid in the chain will be replaced with all 19 other amino acids 

After the addition go to the directory where this repository is and run the line below in terminal to begin the simulation. Make sure that you have installed PyRosetta for that environment. 
```
python3 main.py
```
*Note: make sure that the other lines are commented out. Additionally, note that the simulation can take sometime. During the simulation you will see outputs in your terminal, but you can and should ignore them. A progress bar is also in the output after each iteration is completed and should provide you with an estimated time written on the right hand side of the bar.*

## The following data are the headers of the .`csv` file: 
1. **type**: wild-type or mutant (*only the first row is wild-type*) 
2. **residue_number**: The number of the residue that is being mutated. In "ARVB" sequence, if R is being mutated then the value will be 2 (*one-based indexing*)
3. **previous_aa**: the previous amino acid that was in the wild-type before being mutated at location `residue_number`
4. **new_aa_1l**: the 1 letter code for the new amino acid that will be in the mutant               
5. **new_aa_3l**: the 3 letter code for the new amino acid that will be in the mutant                       
6. **conversion**: Details location and identity of the both former and successor proteins. A String where the value follows the following formula. 
```python
conversion = f"{residue_number}{previous_aa}to{new_aa_1l}"
```
8. **new_seq**: if mutant then the sequence after the mutation otherwise wild-type
9. **fa_score** : Full Atom energy score calculated using `ref2015` formula in PyRosetta              
10. **ddg_score**: the delta delta g value (also calculated using `ref2015`)
11. **hbond_score**: the number of hydrogen-bonds  
12. **sasa_score**: The solvent accessible surface area (SASA) is a measure of the protein's exposure to the solvent. Comparing SASA values between wild-type and mutant proteins can help you determine if the mutation affects the exposure of hydrophobic or hydrophilic residues.            
13. **secondary_structure**: a String for the secondary protein structure
14. **diff_hbonds**: `hbond_score of mutant -  hbond_score of wildtype`
15. **diff_sasa**: `sasa_score of mutant -  sasa_score of wildtype`
16. **diff_secondary_structure**: `secondary_structure of mutant -  secondary_structure of wildtype` calculated by finding out the number of different elements 

