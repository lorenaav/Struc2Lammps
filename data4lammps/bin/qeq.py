"""
Code based on the OpenBabel Qeq implementation available in GitHub
https://github.com/openbabel/openbabel/blob/master/src/charges/qeq.cpp

"""

import numpy as np
from scipy.special import erf
import string

def get_parameters(infile):

    eV = 3.67493245e-2
    Angstrom = 1./0.529177249

    parameters = {}

    with open(infile) as f:
        for i in range(13):
            next(f)
        for line in f:

            data = line.rstrip().split()

            element = data[0]
            elecnegativity = float(data[1])*eV
            hardness = float(data[2])*eV
            sradius = float(data[3])*Angstrom

            basis =  1.0/(sradius*sradius)

            parameters[element] = [elecnegativity, hardness, basis]

    return parameters

def get_elements(infile):

    atomElements = {}

    with open(infile) as f:

        for line in f:

            data = line.rstrip().split()

            typeID = int(data[0])
            typeElement = data[1].split('_')[0]

            atomElements[typeID] = typeElement

    return atomElements

def read_total_charge(inpfile):

    with open(inpfile) as f:

        for line in f:

            if line.startswith('TCHARGE'):

                total_charge = float(line.rstrip().split()[1])

                return total_charge

def calculate_coulomb_integral(a, b, R):
    p = np.sqrt(a*b/(a+b))

    return erf(p*R)/R

def solve_system(C, D):

    charges = np.linalg.solve(C, -D)

    return charges

def fill_J(atoms, J, BasisSet, CoulombMaxDistance):

    nAtoms = len(atoms)

    for k in range(1, nAtoms):
        for l in range(k+1, nAtoms):

            atom1 = atoms[k]
            atom2 = atoms[l]

            xyz_atom1 = np.array(atom1[4:7])
            xyz_atom2 = np.array(atom2[4:7])

            R = np.linalg.norm(xyz_atom1-xyz_atom2)

            if R < CoulombMaxDistance:

                coulomb = calculate_coulomb_integral(BasisSet[k], BasisSet[l], R)

            else:

                coulomb = 1.0/R

            J[k][l] = coulomb
            J[l][k] = coulomb

    for i in range(nAtoms+1):

        J[nAtoms][i] = 1.0

        J[i][nAtoms] = 1.0




def create_C(Q_total, Q_past, fixed, J):


    C = np.zeros((nAtoms, nAtoms))

    for i in range(nAtoms):
        C[0][i] = Q_past[i]

def compute_QEq_charges(atoms, charges_past):

    CoulombThreshold = 1.0e-9

    nAtoms = len(atoms)
    Electronegativity = np.zeros(nAtoms)
    J = np.zeros((nAtoms + 1, nAtoms + 1))
    Voltage = charges_past
    BasisSet = np.zeros(nAtoms)

    parameters = get_parameters('./data4lammps/data/qeq.txt')

    atomElements = get_elements('./types/atom_type.dat')

    for i in range(nAtoms):

        atomType = atoms[i][2]
        atomElement = atomElements[atomType]

        Electronegativity[i] = parameters[atomElement][0]
        J[i][i] = parameters[atomElement][1]
        BasisSet[i] = parameters[atomElement][2]

    total_charge = read_total_charge('./types/atoms.dat')

    print('Total molecule charge: {}'.format(total_charge))

    SmallestGuassianExponent = min(BasisSet)
    CoulombMaxDistance = 2*np.sqrt(-np.log(CoulombThreshold) / SmallestGuassianExponent)

    fill_J(atoms, J, BasisSet, CoulombMaxDistance)

    Voltage[:-1] = Electronegativity
    Voltage[-1] = total_charge

    charges = np.linalg.solve(J, Voltage)

    return charges

def Qeq_charge_equilibration(atoms):

    nAtoms = len(atoms)
    new_charges = np.zeros(nAtoms+1)

    for i in range(1):

        charges_past = new_charges
        new_charges = compute_QEq_charges(atoms, charges_past)

    qeq_charges = new_charges

    return qeq_charges
