import string
import os, sys
import math

#####################################################################################
# Get parameter file name, with installed path
def getPCFFParamFile():

    cwd=os.path.dirname(os.path.abspath(sys.argv[0]))
    datadir=os.path.join(cwd, '..', 'data')

    return os.path.join(datadir, 'PCFF31parameters.txt')
    #return 'PCFF31.txt'

#####################################################################################
# Get atom mass
def getAtommass(atomtypes):

    atommass=[]
    flag=0
    masslist=[]

    Forcefieldfile=getPCFFParamFile()
    fin=open(Forcefieldfile,'r')
    dataline=fin.readline()
    while dataline!="" and dataline!="\n" and flag==0:
        words=dataline[0:len(dataline)-2]
        if str(words).upper()=="ATOMTYPES":
            flag=1
            dataline=fin.readline()
            words=string.split(dataline[0:len(dataline)-1])
            while str(words[0]).upper()!="END":
                atype=str(words[0]).upper()
                amass=eval(words[2])
                for i in range(len(atomtypes)):
                    atomtypeID=atomtypes[i][0]
                    atomtype=atomtypes[i][1]
                    if atype==atomtype.upper():
                        if atomtypeID not in masslist:
                            masslist.append(atomtypeID)
                            atommass.append([atomtypeID,amass,atomtype])
                dataline=fin.readline()
                words=string.split(dataline[0:len(dataline)-1])
        dataline=fin.readline()
    fin.close()

    return atommass
#####################################################################################
# Get bond coeffs
def getBondCoeffs(bondtypes):

    bondcoeffs=[]
    flag=0
    cmark="  #"

    Forcefieldfile=getPCFFParamFile()
    fin=open(Forcefieldfile,'r')
    dataline=fin.readline()
    while dataline!="" and dataline!="\n" and flag==0:
        words=dataline[0:len(dataline)-2]
        if str(words).upper()=="BOND_STRETCH":
            flag=1
            dataline=fin.readline()
            words=string.split(dataline[0:len(dataline)-1])
            while str(words[0]).upper()!="END":
                atype1=str(words[0])
                atype2=str(words[1])

                k2=0.5*eval(words[3])   ##here 1/2*kr*(r-r0)^2, in Lammps using kr*(r-r0)^2
                b0=1.0*eval(words[4])
                if len(words)>5:
                    c3=1.0*eval(words[5])
                    c4=1.0*eval(words[6])
                    k3=k2*c3
                    k4=k2*c4
                else:
                    k3=0.0
                    k4=0.0
                for i in range(len(bondtypes)):
                    bondtypeID=bondtypes[i][0]
                    atom1type=bondtypes[i][1]
                    atom2type=bondtypes[i][2]
                    if atom1type=="c3" or atom1type=="c2":
                        atom1type="c"
                    if atom2type=="c3" or atom2type=="c2":
                        atom2type="c"
                    if (atype1==atom1type and atype2==atom2type) or (atype2==atom1type and atype1==atom2type):
                        bondcoeffs.append([bondtypeID,b0,k2,k3,k4,cmark,atom1type,atom2type])
                dataline=fin.readline()
                words=string.split(dataline[0:len(dataline)-1])
        dataline=fin.readline()
    fin.close()

    # Generate bond coeffs for input scripts
    bondcoeffs_all=[]
    for i in range(len(bondtypes)+1):
        bondcoeffs_all.append([i,0,0,0,0,'#','0','0'])

    for j in range(len(bondcoeffs)):
        bondtypeID=bondcoeffs[j][0]
        bondcoeffs_all[bondtypeID]=bondcoeffs[j]

    return bondcoeffs_all

########################################################################
def outputtypes(title,typelist):

    filename=str(title)+".index"
    fout=open(filename,'w')
    for i in range(len(typelist)):
        outline=" "
        for j in range(len(typelist[i])):
            outline=outline+str(typelist[i][j])+" "
        print >>fout,i+1,outline
    fout.close()

###################################################################
def getAngletypes(angles,atoms,atomtypes):

    angletype=[]
    angletypes=[]
    angletypeID=0

    nangle=len(angles)
    for i in range(nangle):
        latomID=angles[i][2]
        catomID=angles[i][3]
        ratomID=angles[i][4]
        latomtypeID=atoms[latomID-1][2]
        catomtypeID=atoms[catomID-1][2]
        ratomtypeID=atoms[ratomID-1][2]
        if ([latomtypeID,catomtypeID,ratomtypeID] not in angletype) and ([ratomtypeID,catomtypeID,latomtypeID] not in angletype):
            angletypeID=angletypeID+1
            angletype.append([latomtypeID,catomtypeID,ratomtypeID])
            angletypes.append([angletypeID,latomtypeID,catomtypeID,ratomtypeID])
    outputtypes("angletypes",angletypes)
    # update angletype in angles
    for i in range(nangle):
        latomID=angles[i][2]
        catomID=angles[i][3]
        ratomID=angles[i][4]
        latomtypeID=atoms[latomID-1][2]
        catomtypeID=atoms[catomID-1][2]
        ratomtypeID=atoms[ratomID-1][2]
        for j in range(len(angletypes)):
            latomtypeID0=angletypes[j][1]
            catomtypeID0=angletypes[j][2]
            ratomtypeID0=angletypes[j][3]
            if (catomtypeID==catomtypeID0) and ((latomtypeID==latomtypeID0 and ratomtypeID==ratomtypeID0) or (ratomtypeID==latomtypeID0 and latomtypeID==ratomtypeID0)):
                angletypeID=angletypes[j][0]
        angles[i][1]=angletypeID

    for i in range(len(angletypes)):
        latomtypeID=angletypes[i][1]
        catomtypeID=angletypes[i][2]
        ratomtypeID=angletypes[i][3]
        for j in range(len(atomtypes)):
            atomtypeID0=atomtypes[j][0]
            if latomtypeID==atomtypeID0:
                latomtype=atomtypes[j][1]
            if catomtypeID==atomtypeID0:
                catomtype=atomtypes[j][1]
            if ratomtypeID==atomtypeID0:
                ratomtype=atomtypes[j][1]
        angletypes[i][1]=latomtype
        angletypes[i][2]=catomtype
        angletypes[i][3]=ratomtype

    return angletypes,angles

#####################################################################################
# Get angle coeffs
def getAngleCoeffs(angletypes):
 
   anglecoeffs=[]
    flag=0
    cmark="  #"

    Forcefieldfile=getPCFFParamFile()
    fin=open(Forcefieldfile,'r')
    dataline=fin.readline()
    while dataline!="" and dataline!="\n" and flag==0:
        words=dataline[0:len(dataline)-2]
        if str(words).upper()=="ANGLE_BEND":
            flag=1
            dataline=fin.readline()
            words=string.split(dataline[0:len(dataline)-1])
            while str(words[0]).upper()!="END":
                atype1=str(words[0])
                atype2=str(words[1])
                atype3=str(words[2])
                atype=(atype1,atype2,atype3)
                k2=0.5*eval(words[4])
                s0=eval(words[5])
                if len(words)>7:
                    c3=eval(words[6])
                    c4=eval(words[7])
                    k3=k2*c3
                    k4=k2*c4
                else:
                    k3=0
                    k4=0

                for i in range(len(angletypes)):
                    angletypeID=angletypes[i][0]
                    atomtype1=angletypes[i][1]
                    atomtype2=angletypes[i][2]
                    atomtype3=angletypes[i][3]
                    if atomtype1=="c3" or atomtype1=="c2":
                        atomtype1="c"
                    if atomtype2=="c3" or atomtype3=="c2":
                        atomtype2="c"
                    if atomtype3=="c3" or atomtype3=="c2":
                        atomtype3="c"
                    if (atomtype1==atype[0] and atomtype2==atype[1] and atomtype3==atype[2]) or (atomtype3==atype[0] and atomtype2==atype[1] and atomtype1==atype[2]):
                        anglecoeffs.append([angletypeID,s0,k2,k3,k4,cmark,atomtype1,atomtype2,atomtype3])
                dataline=fin.readline()
                words=string.split(dataline[0:len(dataline)-1])
        dataline=fin.readline()
    fin.close()

    # Generate angle coeffs for input scripts
    anglecoeffs_all=[0]
    for i in range(len(angletypes)):
        at1=angletypes[i][1]
        at2=angletypes[i][2]
        at3=angletypes[i][3]        
        anglecoeffs_all.append([i+1,0,0,0,0,'#',at1,at2,at3])

    for j in range(len(anglecoeffs)):
        angletypeID=anglecoeffs[j][0]
        anglecoeffs_all[angletypeID]=anglecoeffs[j]

    return anglecoeffs_all

#####################################################################################
# Get angle coeffs
def getBBCoeffs(angletypes,bondtypes,bondcoeffs):

    BBcoeffs=[0]
    flag=0

    for i in range(len(angletypes)):
        angleID=angletypes[i][0]
        at1=angletypes[i][1]
        at2=angletypes[i][2]
        at3=angletypes[i][3]
        BBcoeffs.append([angleID,0.0,0.0,0.0,"  #",at1,at2,at3])
        for j in range(len(bondtypes)):
            bID=bondtypes[j][0]
            bat1=bondtypes[j][1]
            bat2=bondtypes[j][2]
            if (bat1==at1 and bat2==at2) or (bat2==at1 and bat1==at2):
                b1=bondcoeffs[bID][1]
            if (bat1==at3 and bat2==at2) or (bat2==at3 and bat1==at2):
                b2=bondcoeffs[bID][1]
        BBcoeffs[angleID][2]=b1
        BBcoeffs[angleID][3]=b2

    Forcefieldfile=getPCFFParamFile()
    fin=open(Forcefieldfile,'r')
    dataline=fin.readline()
    while dataline!="" and dataline!="\n" and flag==0:
        words=dataline[0:len(dataline)-2]
        if str(words).upper()=="STRETCH_STRETCH":
            flag=1
            dataline=fin.readline()
            words=string.split(dataline[0:len(dataline)-1])
            while str(words[0]).upper()!="END":
                atype1=str(words[0])
                atype2=str(words[1])
                atype3=str(words[2])
                atype=(atype1,atype2,atype3)
                if len(words)>4:
                    kv=eval(words[4])
                else:
                    kv=0.0
                for i in range(len(angletypes)):
                    angletypeID=angletypes[i][0]
                    atomtype1=angletypes[i][1]
                    atomtype2=angletypes[i][2]
                    atomtype3=angletypes[i][3]
                    if (atomtype1==atype[0] and atomtype2==atype[1] and atomtype3==atype[2]) or (atomtype3==atype[0] and atomtype2==atype[1] and atomtype1==atype[2]):
                        BBcoeffs[angletypeID][1]=kv
                dataline=fin.readline()
                words=string.split(dataline[0:len(dataline)-1])
        dataline=fin.readline()
    fin.close()

    return BBcoeffs

#####################################################################################
# Get bond-angle coeffs
def getBACoeffs(angletypes,bondtypes,bondcoeffs):

    BAcoeffs=[0]
    flag=0

    for i in range(len(angletypes)):
        angleID=angletypes[i][0]
        at1=angletypes[i][1]
        at2=angletypes[i][2]
        at3=angletypes[i][3]
        BAcoeffs.append([angleID,0.0,0.0,0.0,0.0,"  #",at1,at2,at3])
        for j in range(len(bondtypes)):
            bID=bondtypes[j][0]
            bat1=bondtypes[j][1]
            bat2=bondtypes[j][2]
            if (bat1==at1 and bat2==at2) or (bat2==at1 and bat1==at2):
                b1=bondcoeffs[bID][1]
            if (bat1==at3 and bat2==at2) or (bat2==at3 and bat1==at2):
                b2=bondcoeffs[bID][1]
        BAcoeffs[angleID][3]=b1
        BAcoeffs[angleID][4]=b2

    Forcefieldfile=getPCFFParamFile()
    fin=open(Forcefieldfile,'r')
    dataline=fin.readline()
    while dataline!="" and dataline!="\n" and flag==0:
        words=dataline[0:len(dataline)-2]
        if str(words).upper()=="STRETCH_BEND_STRETCH":
            flag=1
            dataline=fin.readline()
            words=string.split(dataline[0:len(dataline)-1])
            while str(words[0]).upper()!="END":
                atype1=str(words[0])
                atype2=str(words[1])
                atype3=str(words[2])
                atype=(atype1,atype2,atype3)
                if len(words)>4:
                    kv1=eval(words[4])
                    kv2=eval(words[5])
                else:
                    kv1=0.0
                    kv2=0.0
                for i in range(len(angletypes)):
                    angletypeID=angletypes[i][0]
                    atomtype1=angletypes[i][1]
                    atomtype2=angletypes[i][2]
                    atomtype3=angletypes[i][3]
                    if (atomtype1==atype[0] and atomtype2==atype[1] and atomtype3==atype[2]):
                        BAcoeffs[angletypeID][1]=kv1
                        BAcoeffs[angletypeID][2]=kv2
                    if (atomtype3==atype[0] and atomtype2==atype[1] and atomtype1==atype[2]):
                        BAcoeffs[angletypeID][1]=kv2
                        BAcoeffs[angletypeID][2]=kv1
                dataline=fin.readline()
                words=string.split(dataline[0:len(dataline)-1])
        dataline=fin.readline()
    fin.close()

    return BAcoeffs

###################################################################
def getDihstypes(dihs,atoms,atomtypes):

    dihstype=[]
    dihstypes=[]
    dihstypeID=0

    ndihs=len(dihs)
    for i in range(ndihs):
        catom1ID=dihs[i][3]
        catom2ID=dihs[i][4]
        eatom1ID=dihs[i][2]
        eatom2ID=dihs[i][5]
        catom1typeID=atoms[catom1ID-1][2]
        catom2typeID=atoms[catom2ID-1][2]
        eatom1typeID=atoms[eatom1ID-1][2]
        eatom2typeID=atoms[eatom2ID-1][2]
        comb1=[eatom1typeID,catom1typeID,catom2typeID,eatom2typeID]
        if (comb1 not in dihstype):
            dihstypeID=dihstypeID+1
            dihstype.append(comb1)
            dihstypes.append([dihstypeID,eatom1typeID,catom1typeID,catom2typeID,eatom2typeID])
    outputtypes("torsiontypes",dihstypes)
    for i in range(len(dihstypes)):
        atom1typeID=dihstypes[i][1]
        atom2typeID=dihstypes[i][2]
        atom3typeID=dihstypes[i][3]
        atom4typeID=dihstypes[i][4]
        for j in range(len(atomtypes)):
            atomtypeID0=atomtypes[j][0]
            if atom1typeID==atomtypeID0:
                atom1type=atomtypes[j][1]
            if atom2typeID==atomtypeID0:
                atom2type=atomtypes[j][1]
            if atom3typeID==atomtypeID0:
                atom3type=atomtypes[j][1]
            if atom4typeID==atomtypeID0:
                atom4type=atomtypes[j][1]

        dihstypes[i][1]=atom1type
        dihstypes[i][2]=atom2type
        dihstypes[i][3]=atom3type
        dihstypes[i][4]=atom4type
    # update dihstype in dihss
    for i in range(ndihs):
        catom1ID=dihs[i][3]
        atom1typeID=atoms[catom1ID-1][2]
        atom1type=atomtypes[atom1typeID-1][1]
        catom2ID=dihs[i][4]
        atom2typeID=atoms[catom2ID-1][2]
        atom2type=atomtypes[atom2typeID-1][1]

        eatom1ID=dihs[i][2]
        atom1typeID=atoms[eatom1ID-1][2]
        atom3type=atomtypes[atom1typeID-1][1]
        eatom2ID=dihs[i][5]
        atom2typeID=atoms[eatom2ID-1][2]
        atom4type=atomtypes[atom2typeID-1][1]

        for j in range(len(dihstypes)):
            atom1type0=dihstypes[j][2]
            atom2type0=dihstypes[j][3]
            atom3type0=dihstypes[j][1]
            atom4type0=dihstypes[j][4]
            if ((atom1type==atom1type0) and (atom2type==atom2type0) and (atom3type==atom3type0) and (atom4type==atom4type0)):
                dihstypeID=dihstypes[j][0]
                break
        dihs[i][1]=dihstypeID

    return dihstypes,dihs

#####################################################################################
# Get dihs coeffs: harmonic
def getDihsCoeffs(dihstypes):

    type_done=[]
    flag=0
    cmark="  #"

    dihscoeffs=[0]
    for i in range(len(dihstypes)):
        dihstypeID=dihstypes[i][0]
        atom1type=dihstypes[i][1]
        atom2type=dihstypes[i][2]
        atom3type=dihstypes[i][3]
        atom4type=dihstypes[i][4]
        dihscoeffs.append([dihstypeID,0.0,0.0,0.0,0.0,0.0,0.0,cmark,atom1type,atom2type,atom3type,atom4type])

    Forcefieldfile=getPCFFParamFile()
    fin=open(Forcefieldfile,'r')
    dataline=fin.readline()
    while dataline!="" and dataline!="\n" and flag==0:
        words=dataline[0:len(dataline)-2]
        if str(words).upper()=="TORSIONS":
            flag=1
            dataline=fin.readline()
            words=string.split(dataline[0:len(dataline)-1])
            while str(words[0]).upper()!="END":
                atype1=str(words[1])
                atype2=str(words[2])
                atype0=str(words[0])
                atype3=str(words[3])

                tp=([0.0,0.0],[0.0,0.0],[0.0,0.0])

                kv=0.5*eval(words[5])
                n=int(eval(words[6]))-1
                cfv=eval(words[7])
                fv=math.acos(cfv)/math.pi*180.0  ##angle
                tp[n][0]=kv
                tp[n][1]=fv
                if len(words)>11:
                    k2=0.5*eval(words[8])
                    n2=int(eval(words[9]))-1
                    cf2=eval(words[10])
                    f2=math.acos(cf2)/math.pi*180.0  ##angle
                    tp[n2][0]=k2
                    tp[n2][1]=f2

                    k3=0.5*eval(words[11])
                    n3=int(eval(words[12]))-1
                    cf3=eval(words[13])
                    f3=math.acos(cf3)/math.pi*180.0  ##angle
                    tp[n3][0]=k3
                    tp[n3][1]=f3
                if len(words)>8 and len(words)<12 and words[8]!='P*':
                    k2=0.5*eval(words[8])
                    n2=int(eval(words[9]))-1
                    cf2=eval(words[10])
                    f2=math.acos(cf2)/math.pi*180.0  ##angle
                    tp[n2][0]=k2
                    tp[n2][1]=f2

                for i in range(len(dihstypes)):
                    dihstypeID=dihstypes[i][0]
                    atom1type=dihstypes[i][1]
                    atom2type=dihstypes[i][2]
                    atom3type=dihstypes[i][3]
                    atom4type=dihstypes[i][4]
                    if ((atype1==atom2type and atype2==atom3type) or (atype1==atom3type and atype2==atom2type)) and ((atype0==atom1type and atype3==atom4type) or (atype0==atom4type and atype3==atom1type)):
                        if dihstypeID not in type_done:
                            type_done.append(dihstypeID)
                            k1=tp[0][0]
                            f1=tp[0][1]
                            k2=tp[1][0]
                            f2=tp[1][1]
                            k3=tp[2][0]
                            f3=tp[2][1]
                            dihscoeffs[dihstypeID][1]=k1
                            dihscoeffs[dihstypeID][2]=f1
                            dihscoeffs[dihstypeID][3]=k2
                            dihscoeffs[dihstypeID][4]=f2
                            dihscoeffs[dihstypeID][5]=k3
                            dihscoeffs[dihstypeID][6]=f3
                dataline=fin.readline()
                words=string.split(dataline[0:len(dataline)-1])
        dataline=fin.readline()
    fin.close()

    return dihscoeffs

#####################################################################################
def getMBTCoeffs(dihstypes,bondtypes,bondcoeffs):

    flag=0
    type_done=[]
    MBTcoeffs=[0]
    for i in range(len(dihstypes)):
        torsionID=dihstypes[i][0]
        at1=dihstypes[i][1]
        at2=dihstypes[i][2]
        at3=dihstypes[i][3]
        at4=dihstypes[i][4]
        MBTcoeffs.append([torsionID,0.0,0.0,0.0,0.0,"  #",at1,at2,at3,at4])
        for j in range(len(bondtypes)):
            bID=bondtypes[j][0]
            bat1=bondtypes[j][1]
            bat2=bondtypes[j][2]
            if (bat1==at2 and bat2==at3) or (bat2==at3 and bat1==at2):
                b1=bondcoeffs[bID][1]
        MBTcoeffs[i+1][4]=b1

    Forcefieldfile=getPCFFParamFile()
    fin=open(Forcefieldfile,'r')
    dataline=fin.readline()
    while dataline!="" and dataline!="\n" and flag==0:
        words=dataline[0:len(dataline)-2]
        if str(words).upper()=="TORSION_STRETCH":
            flag=1
            dataline=fin.readline()
            words=string.split(dataline[0:len(dataline)-1])
            while str(words[0]).upper()!="END":
                atype0=str(words[0])
                atype1=str(words[1])
                atype2=str(words[2])
                atype3=str(words[3])

                if len(words)>5:
                    n1=eval(words[5])
                    A1=eval(words[6])
                    n2=eval(words[7])
                    A2=eval(words[8])
                    n3=eval(words[9])
                    A3=eval(words[10])

                for i in range(len(dihstypes)):
                    dihstypeID=dihstypes[i][0]
                    atom1type=dihstypes[i][1]
                    atom2type=dihstypes[i][2]
                    atom3type=dihstypes[i][3]
                    atom4type=dihstypes[i][4]
                    if ((atype1==atom2type and atype2==atom3type) or (atype1==atom3type and atype2==atom2type)) and ((atype0==atom1type and atype3==atom4type) or (atype0==atom4type and atype3==atom1type)):
                        if dihstypeID not in type_done:
                            type_done.append(dihstypeID)
                            MBTcoeffs[dihstypeID][1]=A1
                            MBTcoeffs[dihstypeID][2]=A2
                            MBTcoeffs[dihstypeID][3]=A3
                dataline=fin.readline()
                words=string.split(dataline[0:len(dataline)-1])
        dataline=fin.readline()
    fin.close()

    return MBTcoeffs

#####################################################################################
def getEBTCoeffs(dihstypes,bondtypes,bondcoeffs):

    flag=0
    type_done=[]
    EBTcoeffs=[0]
    for i in range(len(dihstypes)):
        torsionID=dihstypes[i][0]
        at1=dihstypes[i][1]
        at2=dihstypes[i][2]
        at3=dihstypes[i][3]
        at4=dihstypes[i][4]
        EBTcoeffs.append([torsionID,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,"  #",at1,at2,at3,at4])
        for j in range(len(bondtypes)):
            bID=bondtypes[j][0]
            bat1=bondtypes[j][1]
            bat2=bondtypes[j][2]
            if (bat1==at1 and bat2==at2) or (bat2==at1 and bat1==at2):
                b1=bondcoeffs[bID][1]
            if (bat1==at3 and bat2==at4) or (bat2==at3 and bat1==at4):
                b3=bondcoeffs[bID][1]
        EBTcoeffs[i+1][7]=b1
        EBTcoeffs[i+1][8]=b3

    Forcefieldfile=getPCFFParamFile()
    fin=open(Forcefieldfile,'r')
    dataline=fin.readline()
    while dataline!="" and dataline!="\n" and flag==0:
        words=dataline[0:len(dataline)-2]
        if str(words).upper()=="STRETCH_TORSION_STRETCH":
            flag=1
            dataline=fin.readline()
            words=string.split(dataline[0:len(dataline)-1])
            while str(words[0]).upper()!="END":
                atype0=str(words[0])
                atype1=str(words[1])
                atype2=str(words[2])
                atype3=str(words[3])

                if len(words)>5:
                    n1=eval(words[5])
                    B1=eval(words[6])
                    C1=eval(words[7])
                    n2=eval(words[8])
                    B2=eval(words[9])
                    C2=eval(words[10])
                    n3=eval(words[11])
                    B3=eval(words[12])
                    C3=eval(words[13])

                for i in range(len(dihstypes)):
                    dihstypeID=dihstypes[i][0]
                    atom1type=dihstypes[i][1]
                    atom2type=dihstypes[i][2]
                    atom3type=dihstypes[i][3]
                    atom4type=dihstypes[i][4]
                    if ((atype1==atom2type and atype2==atom3type) or (atype1==atom3type and atype2==atom2type)) and (atype0==atom1type and atype3==atom4type):
                        if dihstypeID not in type_done:
                            type_done.append(dihstypeID)
                            EBTcoeffs[dihstypeID][1]=B1
                            EBTcoeffs[dihstypeID][2]=B2
                            EBTcoeffs[dihstypeID][3]=B3
                            EBTcoeffs[dihstypeID][4]=C1
                            EBTcoeffs[dihstypeID][5]=C2
                            EBTcoeffs[dihstypeID][6]=C3
                    if ((atype1==atom2type and atype2==atom3type) or (atype1==atom3type and atype2==atom2type)) and (atype0==atom4type and atype3==atom1type):
                        if dihstypeID not in type_done:
                            type_done.append(dihstypeID)
                            EBTcoeffs[dihstypeID][1]=C1
                            EBTcoeffs[dihstypeID][2]=C2
                            EBTcoeffs[dihstypeID][3]=C3
                            EBTcoeffs[dihstypeID][4]=B1
                            EBTcoeffs[dihstypeID][5]=B2
                            EBTcoeffs[dihstypeID][6]=B3
                dataline=fin.readline()
                words=string.split(dataline[0:len(dataline)-1])
        dataline=fin.readline()
    fin.close()

    return EBTcoeffs

#####################################################################################
def getATCoeffs(dihstypes,angletypes,anglecoeffs):

    flag=0
    type_done=[]
    ATcoeffs=[0]
    for i in range(len(dihstypes)):
        torsionID=dihstypes[i][0]
        at1=dihstypes[i][1]
        at2=dihstypes[i][2]
        at3=dihstypes[i][3]
        at4=dihstypes[i][4]
        ATcoeffs.append([torsionID,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,"  #",at1,at2,at3,at4])
        for j in range(len(angletypes)):
            bID=angletypes[j][0]
            bat1=angletypes[j][1]
            bat2=angletypes[j][2]
            bat3=angletypes[j][3]
            if bat2==at2 and ((bat3==at1 and bat1==at3) or (bat1==at1 and bat3==at3)):
                f1=anglecoeffs[bID][1]
            if bat2==at3 and ((bat3==at2 and bat1==at4) or (bat1==at2 and bat3==at4)):
                f2=anglecoeffs[bID][1]
        ATcoeffs[i+1][7]=f1
        ATcoeffs[i+1][8]=f2

    Forcefieldfile=getPCFFParamFile()
    fin=open(Forcefieldfile,'r')
    dataline=fin.readline()
    while dataline!="" and dataline!="\n" and flag==0:
        words=dataline[0:len(dataline)-2]
        if str(words).upper()=="BEND_TORSION_BEND":
            flag=1
            dataline=fin.readline()
            words=string.split(dataline[0:len(dataline)-1])
            while str(words[0]).upper()!="END":
                atype0=str(words[0])
                atype1=str(words[1])
                atype2=str(words[2])
                atype3=str(words[3])

                if len(words)>5:
                    n1=eval(words[5])
                    D1=eval(words[6])
                    E1=eval(words[7])
                    n2=eval(words[8])
                    D2=eval(words[9])
                    E2=eval(words[10])
                    n3=eval(words[11])
                    D3=eval(words[12])
                    E3=eval(words[13])

                for i in range(len(dihstypes)):
                    dihstypeID=dihstypes[i][0]
                    atom1type=dihstypes[i][1]
                    atom2type=dihstypes[i][2]
                    atom3type=dihstypes[i][3]
                    atom4type=dihstypes[i][4]
                    if ((atype1==atom2type and atype2==atom3type) or (atype1==atom3type and atype2==atom2type)) and (atype0==atom1type and atype3==atom4type):
                        if dihstypeID not in type_done:
                            type_done.append(dihstypeID)
                            ATcoeffs[dihstypeID][1]=D1
                            ATcoeffs[dihstypeID][2]=D2
                            ATcoeffs[dihstypeID][3]=D3
                            ATcoeffs[dihstypeID][4]=E1
                            ATcoeffs[dihstypeID][5]=E2
                            ATcoeffs[dihstypeID][6]=E3
                    if ((atype1==atom2type and atype2==atom3type) or (atype1==atom3type and atype2==atom2type)) and (atype0==atom4type and atype3==atom1type):
                        if dihstypeID not in type_done:
                            type_done.append(dihstypeID)
                            ATcoeffs[dihstypeID][1]=E1
                            ATcoeffs[dihstypeID][2]=E2
                            ATcoeffs[dihstypeID][3]=E3
                            ATcoeffs[dihstypeID][4]=D1
                            ATcoeffs[dihstypeID][5]=D2
                            ATcoeffs[dihstypeID][6]=D3
                dataline=fin.readline()
                words=string.split(dataline[0:len(dataline)-1])
        dataline=fin.readline()
    fin.close()
    return ATcoeffs

#####################################################################################
def getAATCoeffs(dihstypes,angletypes,anglecoeffs):

    flag=0
    type_done=[]
    AATcoeffs=[0]
    for i in range(len(dihstypes)):
        torsionID=dihstypes[i][0]
        at1=dihstypes[i][1]
        at2=dihstypes[i][2]
        at3=dihstypes[i][3]
        at4=dihstypes[i][4]
        AATcoeffs.append([torsionID,0.0,0.0,0.0,"  #",at1,at2,at3,at4])
        for j in range(len(angletypes)):
            bID=angletypes[j][0]
            bat1=angletypes[j][1]
            bat2=angletypes[j][2]
            bat3=angletypes[j][3]
            if bat2==at2 and ((bat3==at1 and bat1==at3) or (bat1==at1 and bat3==at3)):
                f1=anglecoeffs[bID][1]
            if bat2==at3 and ((bat3==at2 and bat1==at4) or (bat1==at2 and bat3==at4)):
                f2=anglecoeffs[bID][1]
        AATcoeffs[i+1][2]=f1
        AATcoeffs[i+1][3]=f2

    Forcefieldfile=getPCFFParamFile()
    fin=open(Forcefieldfile,'r')
    dataline=fin.readline()
    while dataline!="" and dataline!="\n" and flag==0:
        words=dataline[0:len(dataline)-2]
        if str(words).upper()=="TORSION_BEND_BEND":
            flag=1
            dataline=fin.readline()
            words=string.split(dataline[0:len(dataline)-1])
            while str(words[0]).upper()!="END":
                atype0=str(words[0])
                atype1=str(words[1])
                atype2=str(words[2])
                atype3=str(words[3])

                if len(words)>5:
                    M=eval(words[5])

                for i in range(len(dihstypes)):
                    dihstypeID=dihstypes[i][0]
                    atom1type=dihstypes[i][1]
                    atom2type=dihstypes[i][2]
                    atom3type=dihstypes[i][3]
                    atom4type=dihstypes[i][4]
                    if ((atype1==atom2type and atype2==atom3type) or (atype1==atom3type and atype2==atom2type)) and ((atype0==atom1type and atype3==atom4type) or (atype0==atom4type and atype3==atom1type)):
                        if dihstypeID not in type_done:
                            type_done.append(dihstypeID)
                            AATcoeffs[dihstypeID][1]=M
                dataline=fin.readline()
                words=string.split(dataline[0:len(dataline)-1])
        dataline=fin.readline()
    fin.close()

    return AATcoeffs

#####################################################################################
def getBB13Coeffs(dihstypes,bondtypes,bondcoeffs):

    flag=0
    type_done=[]
    BB13coeffs=[0]
    for i in range(len(dihstypes)):
        torsionID=dihstypes[i][0]
        at1=dihstypes[i][1]
        at2=dihstypes[i][2]
        at3=dihstypes[i][3]
        at4=dihstypes[i][4]
        BB13coeffs.append([torsionID,0.0,0.0,0.0,"  #",at1,at2,at3,at4])
        for j in range(len(bondtypes)):
            bID=bondtypes[j][0]
            bat1=bondtypes[j][1]
            bat2=bondtypes[j][2]
            if (bat1==at1 and bat2==at2) or (bat2==at1 and bat1==at2):
                b1=bondcoeffs[bID][1]
            if (bat1==at3 and bat2==at4) or (bat2==at3 and bat1==at4):
                b3=bondcoeffs[bID][1]
        BB13coeffs[i+1][2]=b1
        BB13coeffs[i+1][3]=b3

    Forcefieldfile=getPCFFParamFile()
    fin=open(Forcefieldfile,'r')
    dataline=fin.readline()
    while dataline!="" and dataline!="\n" and flag==0:
        words=dataline[0:len(dataline)-2]
        if str(words).upper()=="SEPARATED_STRETCH_STRETCH":
            flag=1
            dataline=fin.readline()
            words=string.split(dataline[0:len(dataline)-1])
            while str(words[0]).upper()!="END":
                atype0=str(words[0])
                atype1=str(words[1])
                atype2=str(words[2])
                atype3=str(words[3])

                if len(words)>5:
                    M=eval(words[5])

                for i in range(len(dihstypes)):
                    dihstypeID=dihstypes[i][0]
                    atom1type=dihstypes[i][1]
                    atom2type=dihstypes[i][2]
                    atom3type=dihstypes[i][3]
                    atom4type=dihstypes[i][4]
                    if ((atype1==atom2type and atype2==atom3type) or (atype1==atom3type and atype2==atom2type)) and ((atype0==atom1type and atype3==atom4type) or (atype0==atom4type and atype3==atom1type)):
                        if dihstypeID not in type_done:
                            type_done.append(dihstypeID)
                            BB13coeffs[dihstypeID][1]=M
                dataline=fin.readline()
                words=string.split(dataline[0:len(dataline)-1])
        dataline=fin.readline()
    fin.close()

    return BB13coeffs

###################################################################
def getImpstypes(imps,atoms,atomtypes):

    impstype=[]
    impstypes=[]
    impstypeID=0
    
    nimps=len(imps)
    for i in range(nimps):
        atom1ID=imps[i][2]
        atom2ID=imps[i][3]
        atom3ID=imps[i][4]
        atom4ID=imps[i][5]
        atom1typeID=atoms[atom1ID-1][2]
        atom2typeID=atoms[atom2ID-1][2]
        atom3typeID=atoms[atom3ID-1][2]
        atom4typeID=atoms[atom4ID-1][2]
        if [atom1typeID,atom2typeID,atom3typeID,atom4typeID] not in impstype:
            impstypeID=impstypeID+1
            impstype.append([atom1typeID,atom2typeID,atom3typeID,atom4typeID])
            impstypes.append([impstypeID,atom1typeID,atom2typeID,atom3typeID,atom4typeID])
    # replace atomtypeID with atomtype
    for i in range(len(impstypes)):
        for k in range(1,5):
            atomtypeID=impstypes[i][k]
            for j in range(len(atomtypes)):
                atomtypeID0=atomtypes[j][0]
                if atomtypeID==atomtypeID0:
                    atomtype=atomtypes[j][1]
            impstypes[i][k]=atomtype  
    # update impstype in imps list    
    for i in range(nimps):
        impstype=[]
        for k in range(1,5):
            catomID=imps[i][k]
            atomtypeID=atoms[catomID-1][2]
            atomtype=atomtypes[atomtypeID-1][1]
            impstype.append(atomtype)
        for j in range(len(impstypes)):
            impstype0=[impstypes[j][1],impstypes[j][2],impstypes[j][3],impstypes[j][3]]                
            if impstype==impstype0:
                impstypeID=impstypes[j][0]
        imps[i][1]=impstypeID        
    outputtypes("impropertypes",impstypes)        
 
    return impstypes,imps

#####################################################################################
def getImpsCoeffs(impstypes):

    impscoeffs=[0]
    flag=0

    for i in range(len(impstypes)):
        impID=impstypes[i][0]
        at1=impstypes[i][1]
        at2=impstypes[i][2]
        at3=impstypes[i][3]
        at4=impstypes[i][4]
        impscoeffs.append([impID,0.0,0.0,"  #",at1,at2,at3,at4])

    Forcefieldfile=getPCFFParamFile()
    fin=open(Forcefieldfile,'r')
    dataline=fin.readline()
    while dataline!="" and dataline!="\n" and flag==0:
        words=dataline[0:len(dataline)-2]
        if str(words).upper()=="INVERSIONS":
            flag=1
            dataline=fin.readline()
            words=string.split(dataline[0:len(dataline)-1])
            while str(words[0]).upper()!="END":
                atype1=str(words[0])
                atype2=str(words[1])
                atype3=str(words[2])
                atype4=str(words[3])
                ki=0.5*eval(words[5])
                if len(words)>6:
                    if words[6]!='P*':
                        fi=eval(words[6])
                for i in range(len(impstypes)):
                    impstypeID=impstypes[i][0]
                    atomtype1=impstypes[i][1]
                    atomtype2=impstypes[i][2]
                    atomtype3=impstypes[i][3]
                    atomtype4=impstypes[i][4]
                    if (atype1==atomtype2) and (atomtype1 in (atype2,atype3,atype4)) and (atomtype3 in (atype2,atype3,atype4)) and (atomtype4 in (atype2,atype3,atype4)):
                        impscoeffs[impstypeID][1]=ki
                dataline=fin.readline()
                words=string.split(dataline[0:len(dataline)-1])
        dataline=fin.readline()
    fin.close()

    return impscoeffs

#####################################################################################
def getAACoeffs(impstypes,angletypes,anglecoeffs):

    flag=0
    type_done=[]
    AAcoeffs=[0]
    for i in range(len(impstypes)):
        impID=impstypes[i][0]
        at1=impstypes[i][1]
        at2=impstypes[i][2]
        at3=impstypes[i][3]
        at4=impstypes[i][4]
        AAcoeffs.append([impID,0.0,0.0,0.0,0.0,0.0,0.0,"  #",at1,at2,at3,at4])
        f1 = 0.0
        f2 = 0.0
        f3 = 0.0
        for j in range(len(angletypes)):
            bID=angletypes[j][0]
            bat1=angletypes[j][1]
            bat2=angletypes[j][2]
            bat3=angletypes[j][3]
            if bat2==at2 and ((bat3==at1 and bat1==at3) or (bat1==at1 and bat3==at3)):
                f1=anglecoeffs[bID][1]
            if bat2==at2 and ((bat3==at1 and bat1==at4) or (bat1==at1 and bat3==at4)):
                f2=anglecoeffs[bID][1]
            if bat2==at2 and ((bat3==at3 and bat1==at4) or (bat1==at3 and bat3==at4)):
                f3=anglecoeffs[bID][1]
        AAcoeffs[i+1][4]=f1
        AAcoeffs[i+1][5]=f2
        AAcoeffs[i+1][6]=f3

    Forcefieldfile=getPCFFParamFile()
    fin=open(Forcefieldfile,'r')
    dataline=fin.readline()
    while dataline!="" and dataline!="\n" and flag==0:
        words=dataline[0:len(dataline)-2]
        if str(words).upper()=="BEND_BEND":
            flag=1
            dataline=fin.readline()
            words=string.split(dataline[0:len(dataline)-1])
            while str(words[0]).upper()!="END":
                atype1=str(words[0])
                atype2=str(words[1])
                atype3=str(words[2])
                atype4=str(words[3])
                if len(words)>5:
                    M=eval(words[5])

                for i in range(len(impstypes)):
                    impstypeID=impstypes[i][0]
                    atom1type=impstypes[i][1]
                    atom2type=impstypes[i][2]
                    atom3type=impstypes[i][3]
                    atom4type=impstypes[i][4]
                    if (atype1==atom2type) and (atype2==atom1type and atype3==atom3type and atype4==atom4type):
                        type_done.append(impstypeID)
                        AAcoeffs[impstypeID][1]=M
                    if (atype1==atom2type) and (atype2==atom3type and atype3==atom4type and atype4==atom1type):
                        type_done.append(impstypeID)
                        AAcoeffs[impstypeID][2]=M
                    if (atype1==atom2type) and (atype2==atom4type and atype3==atom1type and atype4==atom3type):
                        type_done.append(impstypeID)
                        AAcoeffs[impstypeID][3]=M
                dataline=fin.readline()
                words=string.split(dataline[0:len(dataline)-1])
        dataline=fin.readline()
    fin.close()
 
    return AAcoeffs

#####################################################################################
# Get bond coeffs
def getPairCoeffs(atomtypes):

    n=len(atomtypes)
    vdws=[]
    flag=0

    Forcefieldfile=getPCFFParamFile()
    fin=open(Forcefieldfile,'r')
    dataline=fin.readline()
    while dataline!="" and dataline!="\n" and flag==0:
        words=dataline[0:len(dataline)-2]
        if str(words).upper()=="DIAGONAL_VDW":
            flag=1
            dataline=fin.readline()
            words=string.split(dataline[0:len(dataline)-1])
            while str(words[0]).upper()!="END":
                atype=str(words[0])
                R0=eval(words[2])
                D0=eval(words[3])
                for i in range(n):
                    atomtypeID=atomtypes[i][0]
                    atomtype=atomtypes[i][1]
                    if atomtype=="c3" or atomtype=="c2":
                        atomtype="c"
                    if (atype==atomtype):
                        vdws.append([atomtypeID,R0,D0,atomtype])
                dataline=fin.readline()
                words=string.split(dataline[0:len(dataline)-1])
        dataline=fin.readline()
    fin.close()
    # Generate pair coeffs for input scripts
    paircoeffs=[]
    R0=[]
    D0=[]
    atype=[]
    for i in range(n+1):
        R0.append(0)
        D0.append(0)
        atype.append(i)

    for j in range(len(vdws)):
        atomtypeID=vdws[j][0]
        R0[atomtypeID]=vdws[j][1]
        D0[atomtypeID]=vdws[j][2]
        atype[atomtypeID]=vdws[j][3]

    # LJ pair
    filename='PCFF_LJpaircoeffs.txt'
    fout=open(filename,'w')
    for i in range(1,n+1):
        epsilon=D0[i]
        sigma=R0[i]
        comment=atype[i]
        print >>fout,"pair_coeff\t",i,epsilon,sigma,"\t  #",comment
    fout.close()

########################################################################
def readPairCoeffs():

    paircoeffs=[0]

    fin=open("PCFF_LJpaircoeffs.txt",'r')
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

####################################################################################
#atomtypes=[[1,'c'],[2,'c_1']]
#bondtypes=[[1,'c','c'],[2,'c_1','o_1']]
#amass=getAtommass(atomtypes)
#getPairCoeffs(atomtypes)
#b=getBondCoeffs(bondtypes)
