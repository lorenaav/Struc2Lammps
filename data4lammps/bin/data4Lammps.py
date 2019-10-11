import os
import string
import handleAtoms
import handleBonds
import Gasteiger
import qeq
import PCFF
from getForcefield import *

def checkAtomtype(inpfile):
    atomtypes=handleAtoms.Atomtypes(inpfile)
    Forcefieldfile=getDreidingParamFile()
    typedefault=('H_','C_3','N_3','O_3','F_','S_3','Cl','I_','Br_')

    atypes=[]
    flag=0

    fin=open(Forcefieldfile,'r')
    dataline=fin.readline()
    while dataline!="" and dataline!="\n" and flag==0:
        words=dataline[0:len(dataline)-2]
        if str(words).upper()=="ATOMTYPES":
            flag=1
            dataline=fin.readline()
            words=string.split(dataline[0:len(dataline)-1])
            while str(words[0]).upper()!="END":
                atype=str(words[0])
                atypes.append(atype)
                dataline=fin.readline()
                words=string.split(dataline[0:len(dataline)-1])
        dataline=fin.readline()
    fin.close()

    anychange="NO"
    changed=[]
    for i in range(len(atomtypes)):
        atomtypeID=atomtypes[i][0]
        atomtype=atomtypes[i][1]
        if atomtype not in atypes:
            anychange="YES"
            for j in range(len(typedefault)):
                deftype=typedefault[j]
                if atomtype[0:2]==deftype[0:2]:
                    atomtypes[i][1]=deftype
                    changed.append([atomtype,deftype])
    for i in range(len(atomtypes)):
        atomtypeID=atomtypes[i][0]
        atomtype=atomtypes[i][1]
        if atomtype not in atypes:
            anychange="YES"
            for j in range(len(typedefault)):
                deftype=typedefault[j]
                if atomtype[0]==deftype[0]:
                    atomtypes[i][1]=deftype
                    changed.append([atomtype,deftype])
    for i in range(len(atomtypes)):
        atomtypeID=atomtypes[i][0]
        atomtype=atomtypes[i][1]
        if atomtype not in atypes:
            anychange="YES"
            deftype="C_3"
            atomtypes[i][1]=deftype
            changed.append([atomtype,deftype])
    if anychange=="YES":
        fout=open("atom_type_reassigned.dat",'w')
        for i in range(len(atomtypes)):
            atomtypeID=atomtypes[i][0]
            atomtype=atomtypes[i][1]
            print >>fout,atomtypeID,atomtype
        fout.close()
    if anychange=="YES":
        wout=open("Datafile_warnings1.txt",'w')
        print >>wout,"## ==== Warning: Force field parameters ===="
        print >>wout,"##  Atom type is re-assigned as following:"
        print >>wout,"##",changed
        print >>wout,"## ==== Warning: Force field parameters ===="
        wout.close()
    return atomtypes,changed

def printCoeffs(fout,ptitle,ptypes,ptypecoeffs):
    print >>fout
    print >>fout,ptitle
    print >>fout
    for i in range(len(ptypes)):
        outline=""
        for j in range(len(ptypecoeffs[i+1])):
            outline=outline+str(ptypecoeffs[i+1][j])+" \t"
        print >>fout,outline

# Write out force field parameters.
def outputDreidingCoeffs(fout,atomtypes,bondtypes,angletypes,dihstypes,impstypes):
    warning1=getPairCoeffs(atomtypes)  # coeffs are in file "paircoeffs.txt"
    paircoeffs=readPairCoeffs()
    print >>fout
    print >>fout,"Pair Coeffs"
    print >>fout
    for i in range(len(paircoeffs)):
        print >>fout,'%3i  %12.6f %12.6f   %s %s' % (paircoeffs[i][0],paircoeffs[i][1],paircoeffs[i][2],paircoeffs[i][3],paircoeffs[i][4])
    print >>fout
    print >>fout,"Bond Coeffs"
    print >>fout
    bondcoeffs,warning2=getBondCoeffs(bondtypes)
    for i in range(len(bondtypes)):
        print >>fout,'%3i  %12.6f  %12.6f  %s%s%s%s' % (bondcoeffs[i][0],bondcoeffs[i][1],bondcoeffs[i][2],str("  # "),bondcoeffs[i][3],str(" "),bondcoeffs[i][4])
    print >>fout
    print >>fout,"Angle Coeffs"
    print >>fout
    anglecoeffs,warning3=getAngleCoeffs(angletypes)
    for i in range(len(angletypes)):
        print >>fout,'%3i  %12.6f  %12.6f  %s%s%s' % (anglecoeffs[i][0],anglecoeffs[i][1],anglecoeffs[i][2],str("  # X "),anglecoeffs[i][3],str(" X "))
    print >>fout
    print >>fout,"Dihedral Coeffs"
    print >>fout
    dihscoeffs,warning4=getDihsCoeffs(dihstypes)
    for i in range(len(dihstypes)):
        print >>fout,'%3i  %12.6f  %3i %3i %s%s%s%s%s%s%s%s' % (dihscoeffs[i][0],dihscoeffs[i][1],dihscoeffs[i][3],dihscoeffs[i][2],str("  # "),dihscoeffs[i][4],str(" "),dihscoeffs[i][5],str(" "),dihscoeffs[i][6],str(" "),dihscoeffs[i][7])
    print >>fout
    print >>fout,"Improper Coeffs"
    print >>fout
    impscoeffs,warning5=getImpsCoeffs(impstypes)
    for i in range(len(impscoeffs)):
        print >>fout,'%3i  %12.6f  %12.6f  %s%s%s' % (impscoeffs[i][0],impscoeffs[i][1],impscoeffs[i][2],str("  #  "),impscoeffs[i][3],str(" X X X "))

    if warning1!="" or warning2!="" or warning3!="" or warning4!="" or warning5!="":
        wout=open("Datafile_warnings2.txt",'w')
        print >>wout,"##",warning1
        print >>wout,"##",warning2
        print >>wout,"##",warning3
        print >>wout,"##",warning4
        print >>wout,"##",warning5
        print >>wout,"## ==== Warning: Force field parameters ====="
        wout.close()

# Write out force field parameters.
def outputPCFFCoeffs(fout,atomtypes,bondtypes,angletypes,dihstypes,impstypes):
    PCFF.getPairCoeffs(atomtypes)  # coeffs are in file "paircoeffs.txt"
    paircoeffs=PCFF.readPairCoeffs()
    printCoeffs(fout,"Pair Coeffs",atomtypes,paircoeffs)

    bondcoeffs=PCFF.getBondCoeffs(bondtypes)
    printCoeffs(fout,"Bond Coeffs",bondtypes,bondcoeffs)

    anglecoeffs=PCFF.getAngleCoeffs(angletypes)
    printCoeffs(fout,"Angle Coeffs",angletypes,anglecoeffs)

    BBcoeffs=PCFF.getBBCoeffs(angletypes,bondtypes,bondcoeffs)
    printCoeffs(fout,"BondBond Coeffs",angletypes,BBcoeffs)

    BAcoeffs=PCFF.getBACoeffs(angletypes,bondtypes,bondcoeffs)
    printCoeffs(fout,"BondAngle Coeffs",angletypes,BAcoeffs)

    dihscoeffs=PCFF.getDihsCoeffs(dihstypes)
    printCoeffs(fout,"Dihedral Coeffs",dihstypes,dihscoeffs)

    MBTcoeffs=PCFF.getMBTCoeffs(dihstypes,bondtypes,bondcoeffs)
    printCoeffs(fout,"MiddleBondTorsion Coeffs",dihstypes,MBTcoeffs)

    EBTcoeffs=PCFF.getEBTCoeffs(dihstypes,bondtypes,bondcoeffs)
    printCoeffs(fout,"EndBondTorsion Coeffs",dihstypes,EBTcoeffs)

    ATcoeffs=PCFF.getATCoeffs(dihstypes,angletypes,anglecoeffs)
    printCoeffs(fout,"AngleTorsion Coeffs",dihstypes,ATcoeffs)

    AATcoeffs=PCFF.getAATCoeffs(dihstypes,angletypes,anglecoeffs)
    printCoeffs(fout,"AngleAngleTorsion Coeffs",dihstypes,AATcoeffs)

    BB13coeffs=PCFF.getBB13Coeffs(dihstypes,bondtypes,bondcoeffs)
    printCoeffs(fout,"BondBond13 Coeffs",dihstypes,BB13coeffs)

    impscoeffs=PCFF.getImpsCoeffs(impstypes)
    printCoeffs(fout,"Improper Coeffs",impstypes,impscoeffs)

    AAcoeffs=PCFF.getAACoeffs(impstypes,angletypes,anglecoeffs)
    printCoeffs(fout,"AngleAngle Coeffs",impstypes,AAcoeffs)

def getFileName():

    fin = open('structure.name','r')
    dataline = fin.readline()
    words = string.split(dataline[0:len(dataline)-1])
    structureName = words[0]

    return structureName

def getForcefield():

    fin = open('forcefield.name','r')
    dataline = fin.readline()
    words = string.split(dataline[0:len(dataline)-1])
    forcefieldName = words[0]

    return forcefieldName

def readPairCoeffs():

    paircoeffs=[]

    fin=open("LJpaircoeffs.txt",'r')
    dataline=fin.readline()
    while dataline!="":
        words=string.split(dataline[0:len(dataline)-1])
        atomtype=eval(words[1])
        D0=eval(words[2])
        R0=eval(words[3])
        C1=words[4]
        C2=words[5]
        paircoeffs.append([atomtype,D0,R0,C1,C2])
        dataline=fin.readline()

    return paircoeffs

# Assumes cubic cell
def createReaxDatafile(forcefield,structureName,xlo, xhi, ylo, yhi, zlo, zhi, chargeMethod):

    inpfile='atom_type.dat'
    if str(forcefield).upper()=="DREIDING":
        atomtypes,anychange=checkAtomtype(inpfile)
    else:
        atomtypes=handleAtoms.Atomtypes(inpfile)
        anychange=[]

    print "Atomtypes total=",len(atomtypes)
    inpfile='atoms.dat'
    baseatoms=handleAtoms.AtomsInfo(inpfile)
    print "Atoms total=",len(baseatoms)
    natomtype=len(atomtypes)
    totalatoms=len(baseatoms)

    atommass=getAtommass(atomtypes)

    ####################################################################
    # Output reaxFF data file for lammps
    datafile=structureName+"_reaxFF.data"
    fout=open(datafile,'w')
    print >>fout,"LAMMPS data file using "+forcefield+" for "+structureName
    print >>fout
    print >>fout,str(totalatoms)+"  atoms"
    print >>fout,str(natomtype)+"  atom types"
    print >>fout
    print >>fout,xlo+" "+xhi+" xlo xhi"
    print >>fout,ylo+" "+yhi+" ylo yhi"
    print >>fout,zlo+" "+zhi+" zlo zhi"
    print >>fout
    print >>fout,"Masses"
    print >>fout
    for i in range(len(atommass)):
        atomtype=atommass[i][2]
        #atomtype=atomtype[0]
        print >>fout,'%3i  %12.6f   %s%s' % (atommass[i][0],atommass[i][1],str("  # "), atomtype)

    ####################################################################
    # Output atom data
    print >>fout
    print >>fout, "Atoms"
    print >>fout
    for i in range(len(baseatoms)):
        dataline=str(baseatoms[i])
        w = string.split(dataline[1:len(dataline)-1],',')
        print >>fout, ('%6d %3d %3d %10.5f %15.8f %15.8f %15.8f' %
                       (eval(w[0]),eval(w[1]),eval(w[2]),eval(w[3]),eval(w[4]),eval(w[5]),eval(w[6])))
    print >>fout
    fout.close()
    print(datafile+" created!")

    return datafile
    ####################################################################

def createDatafile(forcefield,structureName, xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz, chargeMethod):

    inpfile='./types/newatom_type.dat'
    if str(forcefield).upper()=="DREIDING":
        atomtypes,anychange=checkAtomtype(inpfile)
    else:
        atomtypes=handleAtoms.Atomtypes(inpfile)
        anychange=[]
    inpfile='./types/newatoms.dat'
    baseatoms=handleAtoms.AtomsInfo(inpfile)
    # Update bondtype if default types used
    inpfile='./types/newbond_type.dat'
    bondtypes=handleBonds.getBondtypes(inpfile)
    for i in range(len(bondtypes)):
        atom1type=bondtypes[i][1]
        atom2type=bondtypes[i][2]
        for j in range(len(anychange)):
            replaced=anychange[j][0]
            defatype=anychange[j][1]
            if atom1type.upper()==replaced.upper():
                bondtypes[i][1]=defatype
            if atom2type.upper()==replaced.upper():
                bondtypes[i][2]=defatype                
    inpfile='./types/newbonds.dat'
    basebonds=handleBonds.getBonds(inpfile,0,1)
    
    print('Equilibrating charge... \n')

    if chargeMethod == 'Gasteiger':
        # Replace charge for Gasteiger charge
        forcefield=getForcefield()
        Gparas=Gasteiger.getGasteiger_parameters(forcefield)
        Q=Gasteiger.getGasteigerCharge(Gparas,atomtypes,baseatoms,basebonds)

    elif chargeMethod == 'QEq':
        Q = qeq.Qeq_charge_equilibration(baseatoms)

    atomlinks=handleAtoms.AtomLink(baseatoms,basebonds)
    baseangles=handleBonds.createAngles(atomlinks)
    if str(forcefield).upper()=="PCFF":
        angletypes,baseangles=PCFF.getAngletypes(baseangles,baseatoms,atomtypes)
    if str(forcefield).upper()=="DREIDING":
        angletypes,baseangles=handleBonds.getAngletypes(baseangles,baseatoms,atomtypes)

    basedihs=handleBonds.createDihedrals(atomlinks,basebonds)
    if str(forcefield).upper()=="PCFF":
        dihstypes,basedihs=PCFF.getDihstypes(basedihs,baseatoms,atomtypes)
    if str(forcefield).upper()=="DREIDING":
        dihstypes,basedihs=handleBonds.getDihstypes(basedihs,baseatoms,atomtypes)

    baseimps=handleBonds.createImpropers(atomlinks)
    if str(forcefield).upper()=="PCFF":
        impstypes,baseimps=PCFF.getImpstypes(baseimps,baseatoms,atomtypes)
    if str(forcefield).upper()=="DREIDING":
        impstypes,baseimps=handleBonds.getImpstypes(baseimps,baseatoms,atomtypes)

    ####################################################################
    # Total quantities
    natomtype=len(atomtypes)
    nbondtype=len(bondtypes)
    nangletype=len(angletypes)
    ndihstype=len(dihstypes)
    nimpstype=len(impstypes)

    totalatoms=len(baseatoms)
    totalbonds=len(basebonds)
    totalangles=len(baseangles)
    totaldihs=len(basedihs)
    totalimps=len(baseimps)
    ####################################################################
    atommass=getAtommass(atomtypes)
    if str(forcefield).upper()=="PCFF":
        atommass=PCFF.getAtommass(atomtypes)

    # Output Lammps data file
    ####################################################################
    datafile = 'LAMMPSDataFile.data'
    fout=open(datafile,'w')
    print >>fout,"LAMMPS data file using "+forcefield+" for "+structureName
    print >>fout
    print >>fout,str(totalatoms)+"  atoms"
    print >>fout,str(totalbonds)+"  bonds"
    print >>fout,str(totalangles)+"  angles"
    print >>fout,str(totaldihs)+"  dihedrals"
    print >>fout,str(totalimps)+"  impropers"
    print >>fout
    print >>fout,str(natomtype)+"  atom types"
    print >>fout,str(nbondtype)+"  bond types"
    print >>fout,str(nangletype)+"  angle types"
    print >>fout,str(ndihstype)+"  dihedral types"
    print >>fout,str(nimpstype)+"  improper types"
    print >>fout
    print >>fout,xlo+" "+xhi+" xlo xhi"
    print >>fout,ylo+" "+yhi+" ylo yhi"
    print >>fout,zlo+" "+zhi+" zlo zhi"
    if xy == '0.0' and xz == '0.0' and yz == '0.0':
        print >>fout
    else:
        print >>fout,xy+" "+xz+" "+yz+" xy xz yz"
    print >>fout
    print >>fout,"Masses"
    print >>fout
    for i in range(len(atommass)):
        print >>fout,'%3i  %12.6f   %s%s' % (atommass[i][0],atommass[i][1],str("  # "), atommass[i][2])

    ####################################################################
    # Output data
    print >>fout
    print >>fout, "Atoms"
    print >>fout
    for i in range(len(baseatoms)):
        dataline=str(baseatoms[i])
        w = string.split(dataline[1:len(dataline)-1],',')
        print >>fout, ('%6d %3d %3d %10.5f %15.8f %15.8f %15.8f %3d %3d %3d' %
                       (eval(w[0]),eval(w[1]),eval(w[2]),Q[i+1],eval(w[4]),eval(w[5]),eval(w[6]), eval(w[7]), eval(w[8]), eval(w[9])))
    print >>fout
    print >>fout, "Bonds"
    print >>fout
    for i in range(len(basebonds)):
        dataline=str(basebonds[i])
        words = string.split(dataline[1:len(dataline)-1],',')
        outline=""
        for i in range(len(words)):
            outline=outline+str(words[i])+" "
        print >>fout, outline
    print >>fout
    print >>fout, "Angles"
    print >>fout
    for i in range(len(baseangles)):
        dataline=str(baseangles[i])
        words = string.split(dataline[1:len(dataline)-1],',')
        outline=""
        for i in range(len(words)):
            outline=outline+str(words[i])+" "
        print >>fout, outline
    print >>fout
    print >>fout, "Dihedrals"
    print >>fout
    for i in range(len(basedihs)):
        dataline=str(basedihs[i])
        words = string.split(dataline[1:len(dataline)-1],',')
        outline=""
        for i in range(len(words)):
            outline=outline+str(words[i])+" "
        print >>fout, outline
    print >>fout
    print >>fout, "Impropers"
    print >>fout
    for i in range(len(baseimps)):
        dataline=str(baseimps[i])
        words = string.split(dataline[1:len(dataline)-1],',')
        outline=""
        for i in range(len(words)):
            outline=outline+str(words[i])+" "
        print >>fout, outline
    # Coeffs
    if str(forcefield).upper()=="DREIDING":
        outputDreidingCoeffs(fout,atomtypes,bondtypes,angletypes,dihstypes,impstypes)
    if str(forcefield).upper()=="PCFF":
        outputPCFFCoeffs(fout,atomtypes,bondtypes,angletypes,dihstypes,impstypes)
    fout.close()

    print(datafile+" created!")
    if os.path.exists("Datafile_warnings.txt"):
        cmd2="rm Datafile_warnings.txt"
        os.system(cmd2)
    if os.path.exists("Datafile_warnings1.txt"):
        cmd1="cat Datafile_warnings1.txt >>Datafile_warnings.txt"
        cmd2="rm Datafile_warnings1.txt"
        os.system(cmd1)
        os.system(cmd2)
    if os.path.exists("Datafile_warnings2.txt"):
        cmd1="cat Datafile_warnings2.txt >>Datafile_warnings.txt"
        cmd2="rm Datafile_warnings2.txt"
        os.system(cmd1)
        os.system(cmd2)

    return datafile

