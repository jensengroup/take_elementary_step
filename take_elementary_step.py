# written by Jan H. Jensen 
#
from rdkit import Chem
from rdkit.Chem import AllChem
import itertools
from rdkit.Chem import rdmolops
from collections import defaultdict
import copy
import time
import numpy as np
from io import StringIO
import sys
#sio = sys.stderr = StringIO()
from rdkit import rdBase
print(rdBase.rdkitVersion)



def getUA(maxValence_list, valence_list):
    UA = []
    DU = []
    for i, (maxValence,valence) in enumerate(zip(maxValence_list, valence_list)):
        if maxValence - valence > 0:
            UA.append(i)
            DU.append(maxValence - valence)
    return UA,DU


def get_BO(AC,valences):
    BO = AC.copy()
    BO_valence = list(BO.sum(axis=1))
    UA,DU = getUA(valences, BO_valence)

    while len(DU) > 1:
        UA_pairs = list(itertools.combinations(UA, 2))

        for i,j in UA_pairs:
            if BO[i,j] > 0:
                BO[i,j] += 1
                BO[j,i] += 1
                break
        
        BO_valence = list(BO.sum(axis=1))
        UA_new, DU_new = getUA(valences, BO_valence)

        if DU_new != DU:
            UA = copy.copy(UA_new)
            DU = copy.copy(DU_new)
        else:
            break
    
    return BO


def BO_is_OK(BO,AC,charge,DU,atomic_valence_electrons,atomicNumList):

    q = 0
    if heterolytic:
        BO_valences = list(BO.sum(axis=1))
        for i,atom in enumerate(atomicNumList):
            q += get_atomic_charge(atom,atomic_valence_electrons[atom],BO_valences[i])
            if atom == 6:
                number_of_single_bonds_to_C = list(BO[i,:]).count(1)
                if number_of_single_bonds_to_C == 2 and BO_valences[i] == 2:
                    q += 1
                if number_of_single_bonds_to_C == 3 and q + 1 < charge:
                    q += 2
    
    if (BO-AC).sum() == sum(DU) and charge == q:
        return True
    else:
        return False


def get_atomic_charge(atom,atomic_valence_electrons,BO_valence):
    if atom == 1:
        charge = 1 - BO_valence
    elif atom == 5:
        charge = 3 - BO_valence
    elif atom == 15 and BO_valence == 5:
        charge = 0
    elif atom == 16 and BO_valence == 6:
        charge = 0
    else:
        charge = atomic_valence_electrons - 8 + BO_valence
          
    return charge


def clean_charges(mol):
# this is a temporary hack. The real solution is to generate several BO matrices in AC2BO and pick the one
# with the lowest number of atomic charges
#
    rxn_smarts = ['[N+:1]=[*:2]-[O-:3]>>[N+0:1]-[*:2]=[O-0:3]',
                  '[N+:1]=[*:2]-[*:3]=[*:4]-[O-:5]>>[N+0:1]-[*:2]=[*:3]-[*:4]=[O-0:5]']

    for smarts in rxn_smarts:
        patt = Chem.MolFromSmarts(smarts.split(">>")[0])
        while mol.HasSubstructMatch(patt):
            rxn = AllChem.ReactionFromSmarts(smarts)
            ps = rxn.RunReactants((mol,))
            mol = ps[0][0]
                    
    return mol
    


def BO2mol(mol,BO_matrix, atomicNumList,atomic_valence_electrons,mol_charge):
# based on code written by Paolo Toscani

    l = len(BO_matrix)
    l2 = len(atomicNumList)
    BO_valences = list(BO_matrix.sum(axis=1))
    
    if (l != l2):
        raise RuntimeError('sizes of adjMat ({0:d}) and atomicNumList '
            '{1:d} differ'.format(l, l2))

    rwMol = Chem.RWMol(mol)

    bondTypeDict = {
        1: Chem.BondType.SINGLE,
        2: Chem.BondType.DOUBLE,
        3: Chem.BondType.TRIPLE
    }
    
    for i in range(l):
        for j in range(i + 1, l):
            bo = int(round(BO_matrix[i, j]))
            if (bo == 0):
                continue
            bt = bondTypeDict.get(bo, Chem.BondType.SINGLE)
            rwMol.AddBond(i, j, bt)
    mol = rwMol.GetMol()
    
    if heterolytic:
        mol = set_atomic_charges(mol,atomicNumList,atomic_valence_electrons,BO_valences,BO_matrix,mol_charge)
    else:
        mol = set_atomic_radicals(mol,atomicNumList,atomic_valence_electrons,BO_valences)
  
    return mol


def set_atomic_charges(mol,atomicNumList,atomic_valence_electrons,BO_valences,BO_matrix,mol_charge):
    q = 0
    for i,atom in enumerate(atomicNumList):
        a = mol.GetAtomWithIdx(i)
        charge = get_atomic_charge(atom,atomic_valence_electrons[atom],BO_valences[i])
        q += charge
        if atom == 6:
            number_of_single_bonds_to_C = list(BO_matrix[i,:]).count(1)
            if number_of_single_bonds_to_C == 2 and BO_valences[i] == 2:
                    q += 1
                    charge = 0        
            if number_of_single_bonds_to_C == 3 and q + 1 < mol_charge:
                    q += 2
                    charge = 1
 
        if (abs(charge) > 0):
            a.SetFormalCharge(charge)
    #rdmolops.SanitizeMol(mol)
    
    mol = clean_charges(mol)
    
    return mol


def set_atomic_radicals(mol,atomicNumList,atomic_valence_electrons,BO_valences):
# The number of radical electrons = absolute atomic charge
    for i,atom in enumerate(atomicNumList):
        a = mol.GetAtomWithIdx(i)
        charge = get_atomic_charge(atom,atomic_valence_electrons[atom],BO_valences[i])

        if (abs(charge) > 0):
            a.SetNumRadicalElectrons(abs(charge))
            
    return mol


def AC2BO(AC,atomicNumList,charge):
    atomic_valence = defaultdict(list)
    atomic_valence[1] = [1]
    atomic_valence[6] = [4]
    atomic_valence[7] = [4,3]
    atomic_valence[8] = [2,1]
    atomic_valence[9] = [1]
    atomic_valence[14] = [4]
    atomic_valence[15] = [5,4,3]
    atomic_valence[16] = [6,4,2]
    atomic_valence[17] = [1]
    atomic_valence[35] = [1]
    atomic_valence[53] = [1]
    

    atomic_valence_electrons = {}
    atomic_valence_electrons[1] = 1
    atomic_valence_electrons[6] = 4
    atomic_valence_electrons[7] = 5
    atomic_valence_electrons[8] = 6
    atomic_valence_electrons[9] = 7
    atomic_valence_electrons[14] = 4
    atomic_valence_electrons[15] = 5
    atomic_valence_electrons[16] = 6
    atomic_valence_electrons[17] = 7
    atomic_valence_electrons[35] = 7
    atomic_valence_electrons[53] = 7

    valences_list_of_lists = []
    for atomicNum in atomicNumList:
        valences_list_of_lists.append(atomic_valence[atomicNum])

    valences_list = list(itertools.product(*valences_list_of_lists))

    best_BO = AC.copy()

    for valences in valences_list:
        AC_valence = list(AC.sum(axis=1))
        UA,DU_from_AC = getUA(valences, AC_valence)
        if len(UA) == 0 or BO_is_OK(AC,AC,charge,DU_from_AC,atomic_valence_electrons,atomicNumList):
            best_BO = AC.copy()
            break
        else:
            BO = get_BO(AC,valences)
            if BO_is_OK(BO,AC,charge,DU_from_AC,atomic_valence_electrons,atomicNumList):
                best_BO = BO.copy()
                break
            elif BO.sum() > best_BO.sum():
                    best_BO = BO.copy()

    return best_BO,atomic_valence_electrons


def AC2mol(mol,AC,atomicNumList,charge):
    BO,atomic_valence_electrons = AC2BO(AC,atomicNumList,charge)
    mol = BO2mol(mol,BO, atomicNumList,atomic_valence_electrons,charge)
    
    return mol


def get_proto_mol(atomicNumList):
    mol = Chem.MolFromSmarts("[#"+str(atomicNumList[0])+"]")
    rwMol = Chem.RWMol(mol)
    for i in range(1,len(atomicNumList)):
        a = Chem.Atom(atomicNumList[i])
        rwMol.AddAtom(a)
    
    mol = rwMol.GetMol()

    return mol

def get_C(R,num_atoms):
    C_matrices = []
    C_making = []
    C_breaking = []
    shape = (num_atoms,num_atoms)
    for i in xrange(num_atoms):
        for j in xrange(i+1,num_atoms):
            C = np.zeros(shape).astype(int)
            if R[i,j] == 0:
                C[i,j] = 1
                C[j,i] = 1
                C_matrices.append(C)
                C_making.append(C)
            else:
                C[i,j] = -1
                C[j,i] = -1
                C_matrices.append(C)
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
        C_matrices.append(C)
        
    for made in make2:
        C = np.zeros(shape).astype(int)
        for item in made:
            C += C_making[item]
        C_matrices.append(C)
        
    for broken, made in break1_make1:
        C = np.zeros(shape).astype(int)
        C += C_breaking[broken]
        C += C_making[made]
        C_matrices.append(C)
        
    for broken, made in break1_make2:
        C = np.zeros(shape).astype(int)
        C += C_breaking[broken]
        for item in made:
            C += C_making[item]
        C_matrices.append(C)

    for broken, made in break2_make1:
        C = np.zeros(shape).astype(int)
        for item in broken:
            C += C_breaking[item]
        C += C_making[made]
        C_matrices.append(C)
        
    for broken, made in break2_make2:
        C = np.zeros(shape).astype(int)
        for item1 in broken:
            C += C_breaking[item1]
        for item2 in made:
            C += C_making[item2]
        C_matrices.append(C)

    return C_matrices


def get_I_elementary(R,num_atoms):
    I_elementary = []
    C_matrices = get_C(R,num_atoms)
    for C in C_matrices:
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


def take_elementary_step(mol,E_cutoff):
    atomicNumList = [a.GetAtomicNum() for a in mol.GetAtoms()]
    proto_mol = get_proto_mol(atomicNumList)

    AC = Chem.GetAdjacencyMatrix(mol)
    
    num_atoms = len(atomicNumList)
    I_elementary = get_I_elementary(AC,num_atoms)
    
    smiles_list = []
    molecules = []
    raw_smiles_list = []
    raw_molecules = []
    for I in I_elementary:
        if I_is_valid(I,atomicNumList):
            newmol = AC2mol(proto_mol,I,atomicNumList,charge)
            raw_smiles = Chem.MolToSmiles(newmol)
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

    smiles_list.insert(0, Chem.MolToSmiles(mol))
    molecules.insert(0,mol)
    
    return smiles_list, molecules


if __name__ == "__main__":
    smiles_list = ['CC','C=C','C#C']
    smiles_list = ['C=C.C=CC=C']
    smiles_list = ['C(=O)O']

    heterolytic = False
    E_cutoff = 200

    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        rdmolops.Kekulize(mol, clearAromaticFlags = True)
        charge = Chem.GetFormalCharge(mol)
        mol = Chem.AddHs(mol)
        elementary_smiles, elementary_mols = take_elementary_step(mol,E_cutoff)

        print "len(elementary_smiles)",len(elementary_smiles)
        print elementary_smiles




