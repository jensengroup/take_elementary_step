# written by Jan H. Jensen 
#
from rdkit import Chem
from rdkit.Chem import AllChem
import itertools
from rdkit.Chem import rdmolops
from rdkit.Chem import rdFMCS
from collections import defaultdict
import copy
import time
import numpy as np
import xyz2mol
from io import StringIO
import sys
#sio = sys.stderr = StringIO()
#from rdkit import rdBase
#print(rdBase.rdkitVersion)


def get_I_elementary(R,num_atoms,atomicNumList):
    I_elementary = []
    C_making = []
    C_breaking = []
    shape = (num_atoms,num_atoms)
    for i in xrange(num_atoms):
        for j in xrange(i+1,num_atoms):
            C = np.zeros(shape).astype(int)
            if R[i,j] == 0:
                C[i,j] = 1
                C[j,i] = 1
                I = R + C
                if I_is_valid(I,atomicNumList):
                    I_elementary.append(R+C)
                C_making.append(C)
            else:
                C[i,j] = -1
                C[j,i] = -1
                I_elementary.append(R+C)
                C_breaking.append(C)
    
    break1 = range(len(C_breaking))
    make1 = range(len(C_making))
    break2 = list(itertools.combinations(break1,2))
    make2 = list(itertools.combinations(make1,2))

    break1_make1 = list(itertools.product(*[break1,make1]))
    break1_make2 = list(itertools.product(*[break1,make2]))
    break2_make1 = list(itertools.product(*[break2,make1]))
    break2_make2 = list(itertools.product(*[break2,make2]))

    for broken in break2:
        C = np.zeros(shape).astype(int)
        for item in broken:
            C += C_breaking[item]
        I = R + C
        if I_is_valid(I,atomicNumList):
            I_elementary.append(R+C)
        
    for made in make2:
        C = np.zeros(shape).astype(int)
        for item in made:
            C += C_making[item]
        I = R + C
        if I_is_valid(I,atomicNumList):
            I_elementary.append(R+C)
        
    for broken, made in break1_make1:
        C = np.zeros(shape).astype(int)
        C += C_breaking[broken]
        C += C_making[made]
        I = R + C
        if I_is_valid(I,atomicNumList):
            I_elementary.append(R+C)
        
    for broken, made in break1_make2:
        C = np.zeros(shape).astype(int)
        C += C_breaking[broken]
        for item in made:
            C += C_making[item]
        I = R + C
        if I_is_valid(I,atomicNumList):
            I_elementary.append(R+C)

    for broken, made in break2_make1:
        C = np.zeros(shape).astype(int)
        for item in broken:
            C += C_breaking[item]
        C += C_making[made]
        I = R + C
        if I_is_valid(I,atomicNumList):
            I_elementary.append(R+C)
        
    for broken, made in break2_make2:
        C = np.zeros(shape).astype(int)
        for item1 in broken:
            C += C_breaking[item1]
        for item2 in made:
            C += C_making[item2]
        I = R + C
        if I_is_valid(I,atomicNumList):
            I_elementary.append(R+C)


    return I_elementary

def get_BO_energy(m):
    patt = Chem.MolFromSmarts("[*]~[*]")
    bonds = m.GetSubstructMatches(patt)

    bond_types = {}
    bond_types[Chem.BondType.SINGLE] = ""
    bond_types[Chem.BondType.DOUBLE] = "="
    bond_types[Chem.BondType.TRIPLE] = "#"

    bond_energy = {}
    bond_energy["HH"] = 432
    bond_energy["FH"] = 565
    bond_energy["ClH"] = 427
    bond_energy["BrH"] = 363
    bond_energy["HI"] = 295
    
    bond_energy["CH"] = 413
    bond_energy["CC"] = 347
    bond_energy["CN"] = 305
    bond_energy["CO"] = 358
    bond_energy["CF"] = 485
    bond_energy["CCl"] = 339
    bond_energy["BrC"] = 485
    bond_energy["CI"] = 240
    bond_energy["CS"] = 259

    bond_energy["HN"] = 391
    bond_energy["NN"] = 160
    bond_energy["FN"] = 272
    bond_energy["ClN"] = 839
    bond_energy["BrN"] = 243
    bond_energy["NO"] = 201
    bond_energy["HO"] = 467
    bond_energy["OO"] = 146
    bond_energy["FO"] = 190
    bond_energy["ClO"] = 203
    bond_energy["IO"] = 467
    
    bond_energy["FF"] = 154
    bond_energy["ClF"] = 253
    bond_energy["BrF"] = 237
    bond_energy["ClCl"] = 239
    bond_energy["BrCl"] = 218
    bond_energy["BrBr"] = 193
    
    bond_energy["II"] = 149
    bond_energy["ClI"] = 208
    bond_energy["BrI"] = 175
    
    bond_energy["HS"] = 347
    bond_energy["FS"] = 327
    bond_energy["ClS"] = 253
    bond_energy["BrS"] = 218
    bond_energy["SS"] = 266
    
    bond_energy["SiSi"] = 340
    bond_energy["HSi"] = 393
    bond_energy["CSi"] = 360
    bond_energy["OSi"] = 452
    
    bond_energy["C=C"] = 614
    bond_energy["C#C"] = 839
    bond_energy["O=O"] = 495
    bond_energy["C=O"] = 745
    bond_energy["N=O"] = 607
    bond_energy["N=N"] = 418
    bond_energy["N#N"] = 941
    bond_energy["C#N"] = 891
    bond_energy["C=N"] = 615

    energy = 0
    for bond in bonds:
        bond_type = m.GetBondBetweenAtoms(bond[0],bond[1]).GetBondType()
        first_symbol = m.GetAtomWithIdx(bond[0]).GetSymbol() 
        second_symbol = m.GetAtomWithIdx(bond[1]).GetSymbol()
        symbols = [first_symbol,second_symbol]
        symbols.sort()
        patt = symbols[0] + bond_types[bond_type] + symbols[1]
        try:
            energy += bond_energy[patt]
        except:
            energy += 100
    
    return energy


def I_is_valid(I,atomicNumList):
    max_atomic_valence = {}
    max_atomic_valence[1] = 1
    max_atomic_valence[6] = 4
    max_atomic_valence[7] = 4
    max_atomic_valence[8] = 2
    max_atomic_valence[9] = 1
    max_atomic_valence[14] = 4
    max_atomic_valence[15] = 5
    max_atomic_valence[16] = 6
    max_atomic_valence[17] = 1
    max_atomic_valence[35] = 1
    max_atomic_valence[53] = 1
 
    
    valences = list(I.sum(axis=1))
    for atom,valence in zip(atomicNumList,valences):
        if valence > max_atomic_valence[atom]:
            return False
    
    return True

def set_chirality(mol,newmol,atom2chirality):
    Chem.SanitizeMol(newmol)
    chiral = Chem.FindMolChiralCenters(newmol, includeUnassigned=True)
    mol_is_chiral = len(chiral) > 0
    if not mol_is_chiral:
        return newmol

    # Find maximum common substructure (MCS)
    res = rdFMCS.FindMCS([newmol,mol])
    mcs = Chem.MolFromSmarts(res.smartsString)
    mcs_newmol = newmol.GetSubstructMatch(mcs)
    mcs_mol = mol.GetSubstructMatch(mcs)

    atom_map = {key: value for (key, value) in zip(mcs_newmol,mcs_mol)}
  
    for atom, x in chiral:
    
        # Chirality conserved only if entire chiral center (center atom + neighbourgs) matches 
        # that in parent. 
        try:
            neighbourgs = [atom_map[a.GetIdx()] for a in newmol.GetAtomWithIdx(atom).GetNeighbors()]
            neighbourgs_parent = [a.GetIdx() for a in mol.GetAtomWithIdx(atom_map[atom]).GetNeighbors()]
        except:
            continue

        if sorted(neighbourgs) == sorted(neighbourgs_parent):
            if atom2chirality[atom_map[atom]] == 'R':
                newmol.GetAtomWithIdx(atom).SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW)
            else:
                newmol.GetAtomWithIdx(atom).SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW)

    return newmol



def take_elementary_step(mol,charge,E_cutoff,heterolytic,quick):
    chiral_parent = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    parent_is_chiral = len(chiral_parent) > 0
    if parent_is_chiral:
        atom2chirality = {key: value for (key, value) in chiral_parent}

    atomicNumList = [a.GetAtomicNum() for a in mol.GetAtoms()]
    proto_mol = xyz2mol.get_proto_mol(atomicNumList)

    AC = Chem.GetAdjacencyMatrix(mol)
    
    num_atoms = len(atomicNumList)
    I_elementary = get_I_elementary(AC,num_atoms,atomicNumList)
    
    smiles_list = []
    molecules = []
    raw_smiles_list = []
    raw_molecules = []
    for I in I_elementary:
        newmol = xyz2mol.AC2mol(proto_mol,I,atomicNumList,charge,heterolytic,quick)
        if parent_is_chiral:
            newmol = set_chirality(mol,newmol,atom2chirality) 

        raw_smiles = Chem.MolToSmiles(newmol,isomericSmiles=True)
        if raw_smiles not in raw_smiles_list:
            raw_smiles_list.append(raw_smiles)
            raw_molecules.append(newmol)
    
    energy_of_reactant = get_BO_energy(mol)
    for smiles,raw_mol in zip(raw_smiles_list,raw_molecules):
        try:
            test_mol = Chem.MolFromSmiles(smiles)
        except:
            continue
        if test_mol != None:
            energy = get_BO_energy(raw_mol)
            if smiles not in smiles_list and energy_of_reactant-energy < E_cutoff:
                smiles_list.append(smiles)
                molecules.append(raw_mol)

    smiles_list.insert(0, Chem.MolToSmiles(mol,isomericSmiles=True))
    molecules.insert(0,mol)
    
    return smiles_list, molecules


if __name__ == "__main__":
    smiles_list = ['CC','C=C','C#C']
    smiles_list = ['C=C.C=CC=C']
    smiles_list = ['C(=O)O']
    smiles_list = ['C[C@@](C(=C)C=C)(F)O']

    heterolytic = False
    E_cutoff = 200
    quick = True

    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        rdmolops.Kekulize(mol, clearAromaticFlags = True)
        charge = Chem.GetFormalCharge(mol)
        mol = Chem.AddHs(mol)
        elementary_smiles, elementary_mols = take_elementary_step(mol,charge,E_cutoff,heterolytic,quick)

        print "len(elementary_smiles)",len(elementary_smiles)
        print elementary_smiles
