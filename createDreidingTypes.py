import openbabel as ob
import readFiles

def create_dreiding_types(pdb_file, outputname, infile_is_pdb, ffname):

    fftyper = Dreiding()

    path = pdb_file
    s = read_structure(path)
    create_bonds_obabel(s)
    nmol = s.NumResidues()
    natoms = s.NumAtoms()
    nbonds = s.NumBonds()

    #Identify atom types
    
    id_Hbonds = False
    known_atom_types = []
    atom_types = []  # for each atom
    find_atom_types(s, known_atom_types, atom_types, id_Hbonds)

    #Identify bond types
    
    known_bond_types = []
    bond_types = []  # for each bond
    find_bond_types(s, known_bond_types, bond_types, atom_types)

    infile_pdb = './bonds/pdbfile.pdb'
    infile2_lmpdat = './bonds/new.lmpdat'
    infile_lmpdat = './bonds/bonded.lmpdat'
    
    write_atoms_dat(s, atom_types, infile_is_pdb, infile_pdb, infile_lmpdat)
    write_atom_type_dat(known_atom_types, fftyper)
    write_bonds_dat(s, bond_types, infile_pdb, infile2_lmpdat)
    write_bond_type_dat(known_bond_types, known_atom_types, fftyper)

    with open('structure.name', 'w') as f:
        f.write('%s\n' % outputname)

    with open('forcefield.name', 'w') as f:
        f.write('%s\n' % ffname)

def write_bond_type_dat(known_bond_types, known_atom_types, fftyper):

    """Write all bond types to bond_type.dat as input to data4Lammps"""
    with open('./types/bond_type.dat', 'w') as f:
        for i in range(len(known_bond_types)):
            bt = known_bond_types[i]
            ts1 = fftyper.type_str(known_atom_types[bt[0]])
            ts2 = fftyper.type_str(known_atom_types[bt[1]])
            f.write('%d  %s  %s\n' % (i+1, ts1, ts2))

def write_bonds_dat(s, bond_types, infile_pdb, infile2_lmpdat):
    
    # ID_dict1 = create_ID_dictionary1(infile_pdb, infile2_lmpdat)
    """Write all bonds in s to bonds.dat as input to data4Lammps"""
    with open('./types/bonds.dat', 'w') as f, open('./types/bondorder.dat', 'w') as f1:
        f.write("BONDS\n\n")
        for i in range(s.NumBonds()):
            b = s.GetBond(i)
            f.write('%d  %d  %d  %d\n' % (i+1, bond_types[i]+1, \
                                           b.GetBeginAtomIdx(), \
                                           b.GetEndAtomIdx()))
            f1.write('%d  %d\n' % (i+1, b.GetBO()))
            
def write_atom_type_dat(known_atom_types, fftyper):

    """Write all atomic species to atom_type.dat as input to data4Lammps"""
    with open('./types/atom_type.dat', 'w') as f:
        for i in range(len(known_atom_types)):
            at = known_atom_types[i]
            f.write('%d  %s\n' % (i+1, fftyper.type_str(at)))

def write_atoms_dat(s, atom_types, infile_is_pdb, infile_pdb, infile_lmpdat):

    """Write all atoms in s to atoms.dat as input to data4Lammps"""

    images = [0,0,0]
    # ID_dict1 = create_ID_dictionary1(infile_pdb, infile_lmpdat)
    atoms = readFiles.read_atoms_pdb('./bonds/connected_pdb.pdb')

    with open('./types/atoms.dat', 'w') as f, open('./types/atominring.dat','w') as f1:
        natoms = s.NumAtoms()
        f.write('%d atoms\n\nATOMS\n\n' % natoms)

        if not infile_is_pdb:
            id_dict = readFiles.read_atoms_lmpdat('./bonds/new.lmpdat')
            images_dict = readFiles.read_image_info('./bonds/bonded.lmpdat')

        for i in range(natoms):
            a = s.GetAtom(i+1)
            mol = a.GetResidue().GetIdx()+1
            res = a.IsAromatic()
            siz = a.MemberOfRingSize()
            atom_x = a.GetX()
            atom_y = a.GetY()
            atom_z = a.GetZ()

            if not infile_is_pdb:
                images = images_dict[atoms[atom_x, atom_y, atom_z]]
            f.write("{0}  {1}  {2}  {3:.6f}  {4:.6f}  {5:.6f}  {6:.6f} {7} {8} {9}\n".format( \
                 i+1, mol, atom_types[i]+1, 0.0, atom_x, atom_y, atom_z, images[0], images[1], images[2]))           
            if res == True:
                siz = -1    
            f1.write('%d %d\n' % (atoms[atom_x, atom_y, atom_z],siz))

        f.write("\n")
        f.write('TCHARGE {}'.format(s.GetTotalCharge()))

class Typer:
    """Generic forcefield typer"""
    def __init__(self):
        pass

    def type_str(self, at):
        """
        Return a type string for the AtomType at; implemented by derived classes
        """
        pass

class Dreiding(Typer):
    """Apply Dreiding forcefield types to atomic species"""
    def __init__(self):
        Typer.__init__(self)

    def type_str(self, at):
        """Return a Dreiding type string representing the AtomType at"""
        elsym = ob.etab.GetSymbol(at.atomic_num)
        ts = '%s' % elsym
        if len(ts) == 1:
            ts += '_'
        if at.atomic_num == 1 and at.hbond:
            ts += '_HB'
        if elsym not in ['H', 'F', 'Cl', 'Br', 'I', 'Na', 'Ca', 'Fe', 'Zn']:
            # Need hybridization
            if at.is_resonant:
                ts += 'R'
            elif elsym == 'O':
                ts += '%d' % (at.num_bonds+1)
            else:
                ts += '%d' % (at.num_bonds-1)
        return ts

def match_bond_type(bt, t1, t2):
    """Compare bond type tuple bt to (t1,t2) and return True/False"""
    return (bt[0] == t1 and bt[1] == t2) or (bt[0] == t2 and bt[1] == t1)

def find_bond_types(s, known_bond_types, bond_types, atom_types):
    """Identify unique bond types and assign a type to each bond in s"""
    for i in range(s.NumBonds()):
        b = s.GetBond(i)
        t1 = atom_types[b.GetBeginAtomIdx()-1]
        t2 = atom_types[b.GetEndAtomIdx()-1]
        match = False
        ntypes = len(known_bond_types)
        for j in range(ntypes):
            if match_bond_type(known_bond_types[j], t1, t2):
                match = True
                bond_types.append(j)
                break
        if not match:  # new bond type
            bond_types.append(ntypes)
            known_bond_types.append((t1,t2))

class AtomType:
    """Atomic species with hybridization info"""
    def __init__(self, z, nb, res, hb=False):
        self.atomic_num = z
        self.num_bonds = nb
        self.is_resonant = res
        self.hbond = hb

def match_atom_type(at, z, nb, res, hb):
    """Compare AtomType at to atomic info, return True/False"""
    return at.atomic_num == z and at.num_bonds == nb and \
          at.is_resonant == res and at.hbond == hb

def find_atom_types(s, known_atom_types, atom_types, id_Hbonds):
    """Identify unique atomic species and assign a type to each atom in s"""
    for i in range(s.NumAtoms()):
        a = s.GetAtom(i+1)
        z = a.GetAtomicNum()
        nb = a.GetValence()
        res = a.IsAromatic()
        hb = a.IsHbondDonorH() if id_Hbonds else False
        ntypes = len(known_atom_types)
        match = False
        for j in range(ntypes):
            if match_atom_type(known_atom_types[j], z, nb, res, hb):
                match = True
                atom_types.append(j)
                break
        if not match:  # new atom type
            atom_types.append(ntypes)
            known_atom_types.append(AtomType(z, nb, res, hb))

def read_structure(path):

    """Read an input structure file into an OpenBabel molecule"""
    mol = ob.OBMol()
    conv = ob.OBConversion()
    conv.SetInFormat(ob.OBConversion.FormatFromExt(path))
    conv.ReadFile(mol, path)
    p = ob.OBChainsParser()
    p.PerceiveChains(mol)
    return mol

def create_bonds_obabel(s):

    connectivity = connectivity_pdb('./bonds/connected_pdb.pdb')
    atoms = readFiles.read_atoms_pdb('./bonds/connected_pdb.pdb')     
    nAtoms = s.NumAtoms()
    print(nAtoms)
    
    for i in range(nAtoms - 1):
        ai = s.GetAtom(i+1)
        pi = (ai.GetX(), ai.GetY(), ai.GetZ())
        aiID = atoms[pi]
        if not(aiID in connectivity.keys()):
            continue
        if not(aiID in connectivity.keys()):
            continue
        for j in range(i, nAtoms - 1):
            aj = s.GetAtom(j+1)
            pj = (aj.GetX(), aj.GetY(), aj.GetZ())
            ajID = atoms[pj]
            
            if (ajID in connectivity[aiID]):
                s.AddBond(i+1, j+1, 1)
            
            elif (not (ajID in connectivity[aiID])) and (s.GetBond(i+1, j+1) != None):
                ghost_bond = s.GetBond(i+1, j+1)
                s.DeleteBond(ghost_bond)

            elif (not (ajID in connectivity[aiID])) and (s.GetBond(i+1, j+1) != None):
                ghost_bond = s.GetBond(i+1, j+1)
                s.DeleteBond(ghost_bond)


def connectivity_pdb(pdb_file):

    connectivity = {}
    
    with open(pdb_file) as infile:
        for line in infile:
            
             if line.startswith('CONECT'):
                atom_list = list(line.rstrip().split()[1:])
                main_atom = int(atom_list[0])
                bond_list = atom_list[1:]
                bond_list = list(map(int, bond_list))
                connectivity[main_atom] = bond_list

    return connectivity

def find_key(lmpdat_key, pdb_atoms):

    for key in pdb_atoms.keys():

        x_dist = abs(key[0]-lmpdat_key[0])
        y_dist = abs(key[1]-lmpdat_key[1])
        z_dist = abs(key[2]-lmpdat_key[2])        
        difference = max(x_dist, y_dist, z_dist)
        if difference < 0.1:
            return key
            break

def create_ID_dictionary1(infile_pdb, infile_lmpdat):
    
    # ID dict: key = pdb ID
    #          value = lmpdat ID                
    pdb_atoms = readFiles.read_atoms_pdb(infile_pdb)
    lmpdat_atoms = readFiles.read_atoms_lmpdat(infile_lmpdat)    
    ID_dict1 = {}
    
    for key in lmpdat_atoms.keys():
       
        if key not in pdb_atoms.keys():
            new_key = find_key(key, pdb_atoms)
            ID_dict1[pdb_atoms[new_key]] = lmpdat_atoms[key]

        else:
            ID_dict1[pdb_atoms[key]] = lmpdat_atoms[key]

    return ID_dict1
