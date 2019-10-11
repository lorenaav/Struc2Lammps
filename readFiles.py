import os
import subprocess
import openbabel as ob

def convert_structure(infile, opath):

    outfilelmpdat = './bonds/old.lmpdat'
    outfilepdb = './bonds/pdbfile.pdb'
    ovitos = os.path.join(opath, 'ovitos')
    subprocess.call('obabel {0} -O {1}'.format(outfilepdb, './bonds/test.lmpdat'), shell=True)
    ov = subprocess.call('{0} {1}'.format(ovitos, 'getnd.py'), shell=True)
    if ov < 0:
       print('OVITO did not converted file')
    else:
       print('OVITO sucessfully converted file')

    return outfilelmpdat, outfilepdb

def read_n_types(lmpdatFile):

    atomTypes = 0
    with open(lmpdatFile) as f:
        for line in f:
            line = line.lower().rstrip()
            if line.endswith('atom types'):
                atomTypes = line.split()[0]

    return int(atomTypes)

def read_atom_types(lmpdatFile, nAtomTypes):

    types = {}
    typesRead = 0
    with open(lmpdatFile, 'r') as infile:
        for line in infile:
            if line.lower().rstrip() == 'masses':
                next(infile)
                break

        for line in infile:
            typesRead += 1
            if typesRead > nAtomTypes:
                break

            typeID = line.rstrip().split()[0]
            typeName = line.rstrip().split('#')[-1].strip()
            types[typeID] = typeName

    return types

def read_atom_pdb(pdbfile, nAtomTypes):

    types = {}
    ident = []
    with open(pdbfile, 'r') as infile:
        for line in infile:
             if line.startswith('HETATM') or line.startswith('ATOM'):
                atomInfo = line.replace('\n', '')
                atomName = atomInfo[12:16].strip()
                ident.append(atomName)
        
    types1 = []    
    for i in range(len(ident)):
        atomtype = ident[i]
        if atomtype not in types1:
            types1.append(atomtype)   
            types[i+1] = atomtype
    
    return types1, types

def get_mass(types):
    
    m = []
    etab = ob.OBElementTable()
    for value in types:
        a = etab.GetAtomicNum(value)
        c = etab.GetMass(a)
        m.append(c)
        
    return m

def read_bonds_lmpdat(filepath):

    bonds = []
    bondsRead = 0
    with open(filepath, 'r') as infile:
        for line in infile:
            if line.rstrip().endswith(' bonds'):
                nBonds = int(line.split()[0])
            if line.rstrip() == 'Bonds':
                next(infile)
                break
        for line in infile:
            bondsRead += 1
            if bondsRead > nBonds:
                break

            lineElements = line.strip().split()[1:4]
            bondType = eval(lineElements[0])
            atom1 = eval(lineElements[1])
            atom2 = eval(lineElements[2])
            bonds.append([bondType, atom1, atom2])

    return bonds

def read_atoms_pdb(filepath):

    atoms = {}

    with open(filepath, 'r') as infile:
        for line in infile:
            if line.startswith('HETATM') or line.startswith('ATOM'):
                atomInfo = line.replace('\n', '')

                serialNumber = atomInfo[6:11].strip()
                atomName = atomInfo[12:16].strip()
                alternateLocation = atomInfo[16].strip()
                residueName = atomInfo[17:20].strip()
                chainIdentifier = atomInfo[21].strip()
                residueSeqNumber = atomInfo[22:26].strip()
                residueInsertionCode = atomInfo[27].strip()
                xpos = atomInfo[30:38].strip()
                ypos = atomInfo[39:46].strip()
                zpos = atomInfo[47:54].strip()
                occupancy = atomInfo[55:60].strip()
                temperatureFactor = atomInfo[61:66].strip()
                segmentIdentifier = atomInfo[73:76].strip()
                elementSymbol = atomInfo[77:78].strip()
                charge = atomInfo[79:80].strip()

                xpos_float = fix_decimal_places(xpos)
                ypos_float = fix_decimal_places(ypos)
                zpos_float = fix_decimal_places(zpos)

                atoms[(xpos_float, ypos_float, zpos_float)] = int(serialNumber)

    return atoms

def read_atoms_lmpdat(filepath):

    atoms = {}
    bondsRead = 0
    with open(filepath, 'r') as infile:
        for line in infile:
            if line.rstrip().startswith('Atoms'):
                next(infile)
                break

        for line in infile:
            if line == '\n':
                break

            coordinates = line.strip().split()[4:7]
            atomID = line.strip().split()[0]
            coordinates = map(lambda x: fix_decimal_places(x), coordinates)
            atoms[tuple(coordinates)] = int(atomID)

    return atoms

def read_image_info(lmpdat_file):

    images_dict = {}
    bondsRead = 0
    with open(lmpdat_file, 'r') as infile:
        for line in infile:
            if line.rstrip().startswith('Atoms'):
                next(infile)
                break

        for line in infile:
            if line == '\n':
                break

            images = line.strip().split()[7:]
            atomID = line.strip().split()[0]
            images_dict[int(atomID)] = images

    return images_dict

def fix_decimal_places(number):

    fixed_decimals = round(float(number), 3)

    return fixed_decimals
