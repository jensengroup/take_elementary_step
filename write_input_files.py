from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def get_fragments(mol,name):
    fragment_names = []
    fragments = Chem.GetMolFrags(mol,asMols=True)
    labels = ["A","B","C"]
    for label,fragment in zip(labels,fragments):
        fragment_names.append(name+label)
    
    return fragments, fragment_names

def generate_conformations(fragments, max_confs=20):
    for fragment in fragments:
        rot_bond = rdMolDescriptors.CalcNumRotatableBonds(fragment)
        confs = min(3 + 3*rot_bond,max_confs)
        AllChem.EmbedMultipleConfs(fragment,numConfs=confs)
    
    return fragments

def write_mopac_input_file(fragment, fragment_name, keywords):
    number_of_atoms = fragment.GetNumAtoms()
    charge = Chem.GetFormalCharge(fragment)
    symbols = [a.GetSymbol() for a in fragment.GetAtoms()] 
    for i,conf in enumerate(fragment.GetConformers()):
        file_name = fragment_name+"+"+str(i)+".mop"
        with open(file_name, "w") as file:
            file.write(keywords+str(charge)+"\n")
            file.write(" \n")
            file.write(" \n")
            for atom,symbol in enumerate(symbols):
                p = conf.GetAtomPosition(atom)
                line = " ".join((symbol,str(p.x),"1",str(p.y),"1",str(p.z),"1","\n"))
                file.write(line)


def write_xtb_input_file(fragment, fragment_name):
    number_of_atoms = fragment.GetNumAtoms()
    charge = Chem.GetFormalCharge(fragment)
    symbols = [a.GetSymbol() for a in fragment.GetAtoms()] 
    for i,conf in enumerate(fragment.GetConformers()):
        file_name = fragment_name+"+"+str(i)+".xyz"
        with open(file_name, "w") as file:
            file.write(str(number_of_atoms)+"\n")
            file.write("title\n")
            for atom,symbol in enumerate(symbols):
                p = conf.GetAtomPosition(atom)
                line = " ".join((symbol,str(p.x),str(p.y),str(p.z),"\n"))
                file.write(line)
            if charge !=0:
                file.write("$set\n")
                file.write("chrg "+str(charge)+"\n")
                file.write("$end")

# GFN-xTB and MOPAC automatically switches to UHF if the number of electrons is odd, so there is no need
# to specify the multiplicity.
# If you need to do that for another program you can compute the number of electrons by
# atomic_numbers = [a.GetAtomicNum() for a in fragment.GetAtoms()]
# number_of_electrons = sum(atomic_numbers) - charge
        
        

def  write_input_files(mol,name,method,keywords):
# This version writes input files for GFN-xTB and MOPAC. If you want to use another program then replace
# add additional functions
    fragments, fragment_names = get_fragments(mol,name)
    fragments = generate_conformations(fragments)
    for fragment, fragment_name in zip(fragments, fragment_names):
        if method == "xtb":
            write_xtb_input_file(fragment, fragment_name)
        if method == "mopac":
            write_mopac_input_file(fragment, fragment_name,keywords)


if __name__ == "__main__":
    smiles = "CCC.CC"
    name = "test"
    #method = "xtb"
    #keywords = ""
    method = "mopac"
    keywords = "pm3 cycles=200 charge="
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    
    write_input_files(mol,name,method,keywords)
