import sys
import write_input_files as wif
import take_elementary_step as tes
import os
from rdkit import Chem
from rdkit.Chem import rdmolops

# The first entry is the name of the directory that will be created, which also used as the 
# base name of the compounds
# The code currently writes input files for GFN-xTB (basically xyz files). If you want to use another program
# you need to add code in write_input_files.py.
smiles_list = [('diels_alder','C=C.C=CC=C')]
smiles_list = [('formic_acid','C(=O)O')]
#smiles_list = [('ROOR','CC(=O)CO[O].CC(=O)CO[O]')]
#smiles_list = [('ROOR','CC(=O)COOOOCC(=O)C')]

# Make charged fragments, e.g. CH3O- + H+ instead of CH3O + H
charged_fragments = False

# take_elementary_step estimates the molecular energy using bond energies and discards all
# structures whose energy is E_cutoff kcal/mol higher than the starting geometry.
E_cutoff = 100 #kcal/mol

for name,smiles in smiles_list:
    os.mkdir(name)
    os.chdir(name)
    mol = Chem.MolFromSmiles(smiles)
    rdmolops.Kekulize(mol, clearAromaticFlags = True)
    charge = Chem.GetFormalCharge(mol)
    mol = Chem.AddHs(mol)
    elementary_smiles, elementary_mols = tes.take_elementary_step(mol,charge,E_cutoff,charged_fragments)

# The first compound in the lists i the starting struture
    for i,(mol,smiles) in enumerate(zip(elementary_mols,elementary_smiles)):
        new_name = name+str(i)
        wif.write_input_files(mol,new_name)
        print new_name,smiles

