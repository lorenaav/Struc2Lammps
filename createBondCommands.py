import openbabel as ob
import readFiles

def create_bonds_commands(infile, pdbfile, bondscale):

    commands = ''

    nTypes = readFiles.read_n_types(infile)
    nAtomTypes = nTypes
    types, dictt = readFiles.read_atom_pdb(pdbfile, nAtomTypes)
    commands = atom_type_groups_commands(infile, nAtomTypes, types)
    commands += '\n'
    combinations = atom_type_combinations(types)
    cutoff_radii_list = cutoff_radii(combinations, types, bondscale)

    for i in range(len(combinations)):
        type1 = (combinations[i][0])
        type2 = (combinations[i][1])
        cutoff = cutoff_radii_list[i]

        firstLine = 'group    g{0} union {1} {2} \n'.format(i, types[type1], types[type2])
        secondLine = 'fix      f1 g{0} bond/create 1 {1} {2} {3} 1\n'.format(i, type1+1, type2+1, cutoff)
        '''
        if nAngleTypes != 0:
            secondLine += ' atype 1'
        if nDihedralTypes != 0:
            secondLine += ' dtype 1'
        if nImproperTypes != 0:
            secondLine += ' itype 1'
        secondLine += '\n'

        '''

        thirdLine  = 'run      15 \n'
        fourthLine = 'unfix    f1 \n \n'
        commands += (firstLine + secondLine + thirdLine + fourthLine)

    return commands

def atom_type_groups_commands(infile, nTypes, types):
    group_commands = ''
    for i in range(nTypes):
        group_commands += 'group    {0} type {1} \n'.format(types[i], i+1)
    return group_commands

def atom_type_combinations(types):

    atom_type_combinations = []

    for i in range(len(types)):
        for j in range(i, len(types)):
            atom_type_combinations.append((i, j))

    return atom_type_combinations

def cutoff_radii(type_combinations, types, bondscale):

    cutoffs = []
    atomic_numbers = dict(H=1, HE=2, LI=3, BE=4, B=5, C=6, N=7, O=8, F=9, NE=10, NA=11, MG=12,
                          AL=13, SI=14, P=15, S=16, CL=17, AR=18, K=19, CA=20, SC=21, TI=22, V=23, CR=24,
                          MN=25, FE=26, CO=27, NI=28, CU=29, ZN=30, GA=31, GE=32, AS=33, SE=34, BR=35, KR=36,
                   RB=37, SR=38, Y=39, ZR=40, NB=41, MO=42, TC=43, RU=44, RH=45, PD=46, AG=47,
                   CD=48, IN=49, SN=50, SB=51, TE=52, I=53, XE=54, CS=55, BA=56, LA=57,
                   CE=58, PR=59, ND=60, PM=61, SM=62, EU=63, GG=64, TB=65, DY=66, HO=67,
                   ER=68, TM=69, YB=70, LU=71, HF=72, TA=73, W=74, RE=75, OS=76, IR=77,
                   PT=78, AU=79, HG=80, TL=81, PB=82, BI=83, PO=84, AT=85, RN=86, FR=87,
                   RA=88, AC=89, TH=90, PA=91, U=92, NP=93, PU=94, AM=95, CM=96, BK=97,
                   CF=98, ES=99, FM=100, MD=101, NO=102, LR=103, RF=104, DB=105, SG=106, BH=107, HS=108, MT=109)

    for combination in type_combinations:        
        a1Name = types[combination[0]]
        a2Name = types[combination[1]]

        a1Number = atomic_numbers[a1Name]
        a2Number = atomic_numbers[a2Name]

        a1Rad = ob.etab.GetCovalentRad(a1Number)
        a2Rad = ob.etab.GetCovalentRad(a2Number)

        assert type(a1Rad) == float
        assert type(a2Rad) == float
        assert type(bondscale) == float

        cutoffs.append((a1Rad + a2Rad)*bondscale)

    return cutoffs
