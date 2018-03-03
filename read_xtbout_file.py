import numpy as np
import subprocess
import re
import sys
sys.path.insert(0, '/home/jhjensen/github/xyz2mol')
import xyz2mol 

def shell(cmd, shell=False):

    if shell:
        p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    else:
        cmd = cmd.split()
        p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    output, err = p.communicate()
    return output

def get_energy(xtb_out):
    line = shell('grep --text "total E" '+xtb_out+' | tail -1', shell=True)
    energy = re.findall("[-\d]+\.\d+", line)
    if len(energy) != 0:
        energy = energy[0]
        energy = float(energy)
    else:
        energy = 999999.0

    return energy


def read_xtbout_file(xtb_out):
    xyz_coordinates = []
    atomic_symbols = []
    line = shell('grep "charge                     :" '+xtb_out, shell=True)
    charge = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
    charge = int(charge)
    line = shell('grep "number of atoms" '+xtb_out, shell=True)
    number_of_atoms = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
    number_of_atoms = int(number_of_atoms)
# number_of_atoms = 1 exception needs to be implemented
    lines = line = shell('grep -A'+str(number_of_atoms+4)+' "final structure:" '+xtb_out, shell=True)
    lines = lines.split("\n")
    for i in range(4,number_of_atoms+4):
        symbol = lines[i].split(" ")[0].title()
        atomic_symbols.append(symbol)
        xyz = re.findall(r"[-+]?\d*\.\d+|\d+", lines[i])
        xyz = [float(i) for i in xyz]
#        print xyz
        xyz_coordinates.append(xyz)
        
    atomicNumList = xyz2mol.get_atomicNumList(atomic_symbols)
    
    return charge,atomicNumList,xyz_coordinates

if __name__ == "__main__":
    name = "/home/jhjensen/henrik/diels_alder/diels_alder0A+0.xtbout"
    charge,atomicNumList,xyz_coordinates = read_xtbout_file(name)
