import sys
import numpy as np
import subprocess
import openbabel as ob
import readFiles
import createBondCommands
import runData4Lammps

def create_bonds(outfilelmpdat, outfilepdb, bondscale, files_in_lattice, user_boundaries):

    ntypes = readFiles.read_n_types('./bonds/test.lmpdat')
    types, dictt = readFiles.read_atom_pdb(outfilepdb, ntypes)
    mass = readFiles.get_mass(types)
    uboundaries = simulation_box(outfilelmpdat, files_in_lattice, user_boundaries)
    create_lmpdat('./bonds/new.lmpdat', files_in_lattice, uboundaries, outfilepdb)
    write_mass('./bonds/new.lmpdat', ntypes, mass)
    write_coeff('./bonds/new.lmpdat', ntypes)

    create_lammps_in(outfilelmpdat, ntypes, outfilepdb, bondscale)
    return_code = run_lammps()
    if return_code != 0:
        print('LAMMPS did not finish succesfully, can not continue')
        raise Exception('LAMMPS did not finish succesfully, can not continue')
    translate_write_connectivity('./bonds/new.lmpdat', './bonds/bonded.lmpdat', './bonds/pdbfile.pdb')
    print('Bonds have been identified')

def simulation_box(outfilelmpdat, file_with_lattice, user_boundaries):
    
    if file_with_lattice:
        cell_sizes = runData4Lammps.read_cell_sizes(outfilelmpdat)
        print('Box size parameters read:')
        for k, v in cell_sizes.items():
            if len(cell_sizes[k]) < 3:
                print('{0}: {1:.3f} {2:.3f}'.format(k, float(v[0]), float(v[1])))
        uboundaries = user_boundaries
        
    else:
        print('Box size parameters read:')
        print('x: {0:.3f} {1:.3f}'.format(float(user_boundaries[0]), float(user_boundaries[1])))
        print('y: {0:.3f} {1:.3f}'.format(float(user_boundaries[2]), float(user_boundaries[3])))
        print('z: {0:.3f} {1:.3f}'.format(float(user_boundaries[4]), float(user_boundaries[5])))
        
        a = float(user_boundaries[1]); b = float(user_boundaries[3]); c = float(user_boundaries[5])
        alpha = float(user_boundaries[6]); beta = float(user_boundaries[7]); gamma = float(user_boundaries[8])
        x = a
        xy = round(b*(np.cos(np.deg2rad(gamma))),5)
        xz = round(c*(np.cos(np.deg2rad(beta))),5)
        y = np.sqrt(b**2 - xy**2)
        yz = round(b*(c*np.cos(np.deg2rad(alpha))) - xy*xz,5)
        z = np.sqrt(c**2 - xz**2 - yz**2)
        
        boundaries = np.zeros(9)
        boundaries[0] = 0; boundaries[1] = x
        boundaries[2] = 0; boundaries[3] = y
        boundaries[4] = 0; boundaries[5] = z
        boundaries[6] = xy; boundaries[7] = xz; boundaries[8] = yz
        
        return boundaries
  
    
def create_lmpdat(lmpdatFile, file_with_lattice, uboundaries, outfilepdb):

    extra_commands = '''
5 extra bond per atom
10 extra angle per atom
25 extra dihedral per atom
25 extra improper per atom

''' 
    with open('./bonds/old.lmpdat') as oldFile, open('./bonds/new.lmpdat', 'w') as newFile:

        newFile.write('LAMMPS data file \n')
        newFile.write('\n')

        for line in oldFile:
            if line.endswith(' atoms\n'):
                newFile.write(line)
            elif line.endswith(' atom types\n'):
                newFile.write(line)
                newFile.write('\n0 bonds\n1 bond types\n')              
                newFile.write(extra_commands)
                
            elif line.endswith('xhi\n'):
                if file_with_lattice:
                    newFile.write(line)
                else:
                    newFile.write('{0} {1} xlo xhi\n'.format(uboundaries[0], uboundaries[1]))
            elif line.endswith('yhi\n'):
                if file_with_lattice:
                    newFile.write(line)
                else:
                    newFile.write('{0} {1} ylo yhi\n'.format(uboundaries[2], uboundaries[3]))
            elif line.endswith('zhi\n'):
                if file_with_lattice:
                    newFile.write(line)
                    break
                else:
                    newFile.write('{0} {1} zlo zhi\n\n'.format(uboundaries[4], uboundaries[5]))
                    if uboundaries[6] != 0.0 or uboundaries[7] != 0.0 or uboundaries[8] != 0.0:
                        newFile.write('{0} {1} {2} xy xz yz \n\n'.format(uboundaries[6], uboundaries[7], uboundaries[8] ))
                    break
                    
        for line in oldFile:
            print(line)
            if line == '\n':
                break
            else:
                newFile.write(line)

        for line in oldFile:
            newFile.write('\n')
            if line.startswith('Atoms'):
                newFile.write(line)
                newFile.write('\n')
                next(oldFile)
                break
            else:
                newFile.write(line)

        for line in oldFile:
            if line == '\n':
                break
            else:
                newFile.write(line)

def write_mass(lmpdatFile,ntypes, mass):
    with open(lmpdatFile, 'a') as f:
        f.write('\nMasses \n \n')
        for i in range(ntypes):
            f.write('{0} {1} \n'.format(i+1, mass[i]))
            
def write_coeff(lmpdatFile, ntypes):
    with open(lmpdatFile, 'a') as f:
        f.write('\nPair Coeffs \n \n')
        for i in range(ntypes):
            f.write('{}    0.100000     2.00000 \n'.format(i+1))

        f.write('\nBond Coeffs \n \n')
        f.write('1    350.000     1.00000 \n')

def create_lammps_in(lmpdatFile, ntypes, outfilepdb, bondscale):

    constant_string1 = '''###### LAMMPS Fix/Bond Create input file #######
units          real
atom_style     full
boundary       p p p
pair_style     lj/smooth 10.0 12.0 #dummy pair style
special_bonds  lj 0.0 1.0 1.0 extra 10000
bond_style     harmonic
angle_style    harmonic
dihedral_style harmonic
improper_style harmonic

box            tilt large
read_data      ./bonds/new.lmpdat
timestep       1
thermo         1
thermo_style   custom step etotal ke pe temp density vol pxx pyy pzz lx ly lz

'''

    atomCombinations = create_atom_combinations(ntypes)
    outInfile ='./bonds/bondcreate.in'
    with open(outInfile, 'w') as f:
        f.write(constant_string1)
        f.write(createBondCommands.create_bonds_commands(lmpdatFile, outfilepdb, bondscale))
        f.write('write_data ./bonds/bonded.lmpdat \n')

def create_atom_combinations(nAtoms):

    combinations = []
    for i in range(nAtoms):
        for j in range(i, nAtoms):
            combinations.append((i+1,j+1))

    return combinations

def run_lammps():

    return_code = subprocess.Popen('lmp_serial < ./bonds/bondcreate.in > out', shell=True, 
                                   stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
    stdout,stderr = return_code.communicate()

    return return_code

def create_connectivity(infile):

    links_dict = {}
    links = []   
    atoms = readFiles.read_atoms_lmpdat(infile)
    bonds = readFiles.read_bonds_lmpdat(infile)
    print("Number of atoms: ", len(atoms))
    print("Number of bonds: ",len(bonds))
    
    for i in range(len(bonds)):
        atom1 = bonds[i][1]
        atom2 = bonds[i][2]
        links_dict.setdefault(atom1, []).append(atom2)
        links_dict.setdefault(atom2, []).append(atom1)

    for key in sorted(links_dict.keys()):
        tmpLinks = links_dict[key]
        links.append(tmpLinks)

    if len(links_dict) < len(atoms):
        print(len(atoms), len(links_dict))
        print('WARNING: Some atoms are not connected, increasing the bondscale might change this. Data file will be still generated')
    
    return links_dict   

def translate_write_connectivity(infile_lmpdat, connected_infile_lmpdat, infile_pdb):

    ID_dict = create_ID_dictionary(infile_pdb, infile_lmpdat)
    connectivity_lmpdat = create_connectivity(connected_infile_lmpdat)
    connectivity_pdb = {}

    for lmpdat_key, lmpdat_connections in connectivity_lmpdat.items():
        pdb_key = ID_dict[lmpdat_key]
        pdb_connections = map(lambda x: ID_dict[x], lmpdat_connections)
        connectivity_pdb[pdb_key] = pdb_connections

    with open('./bonds/pdbfile.pdb') as oldFile, open('./bonds/connected_pdb.pdb', 'w') as newFile:
        for line in oldFile:
            if line.startswith('ATOM') or line.startswith('HETATM') or line.startswith('COMPND'):
                newFile.write(line)
            if line.startswith('MASTER'):
                masterLine = line

        for key in sorted(connectivity_pdb.keys()):
            connected_atoms = map(str, connectivity_pdb[key])
            connected_atoms = map(lambda x: '{:>5}'.format(x), connected_atoms)
            newFile.write('CONECT   {0:>5} {1}\n'.format(key, ''.join(connected_atoms)))
        newFile.write(masterLine + 'END')

def find_key(pdb_key, lmpdat_atoms):

    for key in lmpdat_atoms.keys():
        x_dist = abs(key[0]-pdb_key[0])
        y_dist = abs(key[1]-pdb_key[1])
        z_dist = abs(key[2]-pdb_key[2])
        difference = max(x_dist, y_dist, z_dist)

        if difference < 0.1:
            return key
            break

def create_ID_dictionary(infile_pdb, infile_lmpdat):

    # ID dict: key = lmpdat ID
    #          value = pdb ID
    pdb_atoms = readFiles.read_atoms_pdb(infile_pdb)
    lmpdat_atoms = readFiles.read_atoms_lmpdat(infile_lmpdat)
    ID_dict = {}

    for key in pdb_atoms.keys():
        if key not in lmpdat_atoms.keys():
            new_key = find_key(key, lmpdat_atoms)
            if new_key == None:
                sys.stderr.write('Coordinates in pdbfile.pdb and bonded.lmpdat are different')
            ID_dict[lmpdat_atoms[new_key]] = pdb_atoms[key]
        else:
            ID_dict[lmpdat_atoms[key]] = pdb_atoms[key]

    return ID_dict
