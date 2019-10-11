import os
import os.path
import shutil
import numpy as np
import subprocess
import readFiles
import createBonds
import createDreidingTypes
import runData4Lammps

# Bonding information in pdb file
def conect_in_pdb(infile):
    
    with open(infile) as pdb_file:    
        for line in pdb_file: 
            if line.startswith('CONECT'):        
                return True
            
        return False

# Initial structure_convert
def structure(name, extension, outfilePDB):

    if extension == '.cif':
       subprocess.call('obabel {0} -O {1} --fillUC strict --title convert ---errorlevel 1'.format(name, outfilePDB), shell=True)
    else:
       subprocess.call('obabel {0} -O {1} --title convert ---errorlevel 1'.format(name, outfilePDB), shell=True)

# Simulation box parameters
def box_params(box_parameters):

    if box_parameters == True:

       boundaries = np.zeros(9)
       return True, boundaries

    else:

       boundaries = np.zeros(9)
       b_file = open('box.txt', 'r')
       cell   = b_file.readline().split()
       angles = b_file.readline().split()
       print(cell, angles)

       boundaries[0] = float(cell[0])
       boundaries[1] = float(cell[1])
       boundaries[2] = float(cell[2])
       boundaries[3] = float(cell[3])
       boundaries[4] = float(cell[4])
       boundaries[5] = float(cell[5])
       boundaries[6] = float(angles[0])
       boundaries[7] = float(angles[1])
       boundaries[8] = float(angles[2])

       return False, boundaries

def run_all(infile, box_parameters, bondscale, ffield, charge, opath):
    
    lattice_in_file, boundaries = box_params(box_parameters)
    extension = os.path.splitext(infile)[1]
    print(lattice_in_file, boundaries, extension)
    extension = os.path.splitext(infile)[1]
    infile_is_pdb = False
    structure(infile, extension, './bonds/pdbfile.pdb') 
 
    if extension == '.pdb' and conect_in_pdb(infile): 
        print('Connectivity already in {}. Copying information ...'.format(infile))
        shutil.copy(infile, './bonds/connected_pdb.pdb')
        outfilelmpdat, outfilepdb = readFiles.convert_structure(infile, opath)
        boundaries = createBonds.simulation_box(outfilelmpdat, lattice_in_file, boundaries)
        infile_is_pdb = True
        
    else:
        print('Creating bonds ...')
        outfilelmpdat, outfilepdb = readFiles.convert_structure(infile, opath)
        createBonds.create_bonds(outfilelmpdat, outfilepdb, bondscale, lattice_in_file, boundaries)

    outputName = infile.split('/')[-1].split('.')[0]    
    print(outputName)    
    createDreidingTypes.create_dreiding_types('./bonds/connected_pdb.pdb', outputName, infile_is_pdb, ffield)
    print('Getting information...')
    
    new_typing_command, data4lammps_command = runData4Lammps.run_data4lammps(charge, infile_is_pdb, boundaries, ffield)

    subprocess.call(new_typing_command)
    subprocess.call(data4lammps_command)
    subprocess.call('rm *index *name')


file = open('input.args', 'r')
name           = str(file.readline().split()[0])
box_parameters = eval(file.readline().split()[0])
bondscale      = float(file.readline().split()[0])
ffield         = str(file.readline().split()[0])
charge         = str(file.readline().split()[0])

print('User defined parameters:')
print('File name: ', name)
print('Box parameters in file?: ', box_parameters)
print('Input values (bondscale, force field, charge eq.): ', bondscale, ffield, charge)

opath = raw_input("Enter the path of your ovitos interpreter: ")
assert os.path.exists(opath), "I did not find the package at, "+str(opath)

run_all(name, box_parameters, bondscale, ffield, charge, opath)
