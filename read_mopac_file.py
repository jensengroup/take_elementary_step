import subprocess
import re
import sys
import xyz2mol 


def shell(cmd, shell=False):

    if shell:
        p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    else:
        cmd = cmd.split()
        p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    output, err = p.communicate()
    return output

def get_energy(mopac_out):
    au_to_kcal = 627.51
    line = shell('grep --text "HEAT OF FORMATION" '+mopac_out+' | tail -1', shell=True)
    energy = re.findall("[-\d]+\.\d+", line)
    if len(energy) != 0:
        energy = energy[0]
        energy = float(energy)
        energy = energy/au_to_kcal 
    else:
        line = shell('grep --text "HEAT" '+mopac_out+' | tail -1', shell=True)
        energy = re.findall("[-\d]+\.\d+", line)
        if len(energy) != 0:
            energy = energy[-1]
            energy = float(energy)
            energy = energy/au_to_kcal
        else: 
            energy = 999999.0

    return energy


def read_mopac_file(mopac_out):
    xyz_coordinates = []
    atomic_symbols = []
    line = shell('grep "CHARGE ON SYSTEM" '+mopac_out, shell=True)
    charge = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
    charge = int(charge)
    line = shell('grep "Empirical Formula:" '+mopac_out, shell=True)
    number_of_atoms = re.findall(r"[-+]?\d*\.\d+|\d+", line)[-1]
    number_of_atoms = int(number_of_atoms)

    special_case_test = shell('grep "CURRENT VALUE OF GRADIENT NORM" '+mopac_out, shell=True)
    special_case = len(special_case_test) > 0
    if special_case:
        lines = shell('grep -A'+str(number_of_atoms+7)+' "CURRENT VALUE OF GEOMETRY" '+mopac_out+' | tail -'+str(number_of_atoms+7), shell=True)
        lines = lines.split("\n")
        for i in range(4,number_of_atoms+4):
            symbol = re.findall(r'\b\w+\b', lines[i])[0]
            atomic_symbols.append(symbol)
            xyz = re.findall(r"[-+]?\d*\.\d+|\d+", lines[i])
            xyz.pop(1)
            xyz.pop(2)
            xyz.pop(3)
            xyz = [float(j) for j in xyz]
            #print xyz
            xyz_coordinates.append(xyz)
    else:
        lines = shell('grep -A'+str(number_of_atoms+1)+' "CARTESIAN COORDINATES" '+mopac_out+' | tail -'+str(number_of_atoms), shell=True)
        lines = lines.split("\n")
        for i in range(number_of_atoms):
            #print re.findall(r"[A-Za-z]", lines[i])
            symbol = re.findall(r'\b\w+\b', lines[i])[1]
            atomic_symbols.append(symbol)
            xyz = re.findall(r"[-+]?\d*\.\d+|\d+", lines[i])
            xyz.pop(0)
            xyz = [float(j) for j in xyz]
            #print xyz
            xyz_coordinates.append(xyz)
           
    atomicNumList = xyz2mol.get_atomicNumList(atomic_symbols)
        
    return charge,atomicNumList,xyz_coordinates


if __name__ == "__main__":
    #name = "finalheat.out"
    #name = "currentbest.out"
    name = "currentvalue.out"
    energy = get_energy(name)
    #print energy
    charge,atomicNumList,xyz_coordinates = read_mopac_file(name)

