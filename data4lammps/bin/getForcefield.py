import os 
import sys
import string
from math import exp

# Get parameter file name, with installed path
def getDreidingParamFile():
    cwd=os.path.dirname(os.path.abspath(sys.argv[0]))
    datadir=os.path.join(cwd, '..', 'data')

    return os.path.join(datadir, 'DreidingX6parameters.txt')
    #return 'DreidingX6parameters.txt'

# Get atom mass
def getAtommass(atomtypes):

    atommass=[]
    flag=0

    Forcefieldfile=getDreidingParamFile()
    fin=open(Forcefieldfile,'r')
    dataline=fin.readline()
    while dataline!="" and dataline!="\n" and flag==0:
        words=dataline[0:len(dataline)-2]
        if str(words).upper()=="ATOMTYPES":
            flag=1
            dataline=fin.readline()
            words=str.split(dataline[0:len(dataline)-1])
            while str(words[0]).upper()!="END":
                atype=str(words[0])
                amass=eval(words[2])
                for i in range(len(atomtypes)):
                    atomtypeID=atomtypes[i][0]
                    atomtype=atomtypes[i][1]
                    if atype==atomtype:
                        atommass.append([atomtypeID,amass,atomtype])
                dataline=fin.readline()
                words=str.split(dataline[0:len(dataline)-1])
        dataline=fin.readline()
    fin.close()
    # if not assigned
    assigned=[]
    for j in range(len(atommass)):
        atomtypeID=atommass[j][0]
        assigned.append(atomtypeID)
    for i in range(len(atomtypes)):
        atomtypeID=atomtypes[i][0]
        atomtype=atomtypes[i][1]
        if atomtypeID not in assigned:
            atype=atomtype[0:2]
            for j in range(len(atommass)):
                btypeID=atommass[j][0]
                amass=atommass[j][1]
                btype=atommass[j][2]
                if (atype==btype[0:2]):
                    atommass.append([atomtypeID,amass,atomtype]) 
                    break
 
    return atommass

# Get bond coeffs
def getBondCoeffs(bondtypes):

    warning=""
    bondcoeffs=[]
    flag=0

    Forcefieldfile=getDreidingParamFile()
    fin=open(Forcefieldfile,'r')
    dataline=fin.readline()
    while dataline!="" and dataline!="\n" and flag==0:
        words=dataline[0:len(dataline)-2]
        if str(words).upper()=="BOND_STRETCH":
            flag=1
            dataline=fin.readline()
            words=str.split(dataline[0:len(dataline)-1])
            while str(words[0]).upper()!="END":
                atype1=str(words[0])
                atype2=str(words[1])

                kr=0.5*eval(words[3])   ##here 1/2*kr*(r-r0)^2, in Lammps using kr*(r-r0)^2
                r0=eval(words[4])

                for i in range(len(bondtypes)):
                    bondtypeID=bondtypes[i][0]
                    atom1type=bondtypes[i][1]
                    atom2type=bondtypes[i][2]
                    if atom1type=="H__HB":
                        atom1type="H___b"
                    if atom2type=="H__HB":
                        atom2type="H___b"
                    if atom1type=="H__HA":
                        atom1type="H___A"
                    if atom2type=="H__HA":
                        atom2type="H___A"
                    if (atype1==atom1type and atype2==atom2type) or (atype2==atom1type and atype1==atom2type):
                        bondcoeffs.append([bondtypeID,kr,r0,atom1type,atom2type])
                dataline=fin.readline()
                words=str.split(dataline[0:len(dataline)-1])
        dataline=fin.readline()
    fin.close()
    # if not assigned
    assigned=[]
    for j in range(len(bondcoeffs)):
        bondtypeID=bondcoeffs[j][0]
        assigned.append(bondtypeID)
    for i in range(len(bondtypes)):
        bondtypeID=bondtypes[i][0]
        atom1type=bondtypes[i][1]
        atom2type=bondtypes[i][2]
        if bondtypeID not in assigned:
            warning=warning+" Bondtype "+str(bondtypeID)+" might be not assigned correctly."
            flag=0
            a1type=atom1type[0]
            a2type=atom2type[0]
            for j in range(len(bondcoeffs)):
                btypeID=bondcoeffs[j][0]
                kr=bondcoeffs[j][1]
                r0=bondcoeffs[j][2]
                b1type=bondcoeffs[j][3]
                b2type=bondcoeffs[j][4]
                if (a1type==b1type[0] and a2type==b2type[0]) or (a1type==b2type[0] and a2type==b1type[0]):
                    bondcoeffs.append([bondtypeID,kr,r0,atom1type,atom2type]) 
                    flag=1
                    break
            if flag==0:
                bondcoeffs.append([bondtypeID,350.0,1.40,atom1type,atom2type])  
    # sorting    
    bondcoeffs.sort()

    return bondcoeffs,warning

# Get bond coeffs
def getAngleCoeffs(angletypes):

    warning=""
    anglecoeffs=[]
    flag=0

    Forcefieldfile=getDreidingParamFile()
    fin=open(Forcefieldfile,'r')
    dataline=fin.readline()
    while dataline!="" and dataline!="\n" and flag==0:
        words=dataline[0:len(dataline)-2]
        if str(words).upper()=="ANGLE_BEND":
            flag=1
            dataline=fin.readline()
            words=str.split(dataline[0:len(dataline)-1])
            while str(words[0]).upper()!="END":
                atype=str(words[1])
                ksita=0.5*eval(words[4])
                sita0=eval(words[5])

                for i in range(len(angletypes)):
                    angletypeID=angletypes[i][0]
                    atomtype=angletypes[i][1]
                    if (atype==atomtype):
                        anglecoeffs.append([angletypeID,ksita,sita0,atomtype])
                dataline=fin.readline()
                words=str.split(dataline[0:len(dataline)-1])
        dataline=fin.readline()
    fin.close()
    # if not assigned
    assigned=[]
    for j in range(len(anglecoeffs)):
        angletypeID=anglecoeffs[j][0]
        assigned.append(angletypeID)
    for i in range(len(angletypes)):
        angletypeID=angletypes[i][0]
        atomtype=angletypes[i][1]
        if angletypeID not in assigned:
            warning=warning+" Angletype "+str(angletypeID)+" might be not assigned correctly."
            flag=0
            atype=atomtype[0]
            for j in range(len(anglecoeffs)):
                btypeID=anglecoeffs[j][0]
                ksita=anglecoeffs[j][1]
                sita0=anglecoeffs[j][2]
                btype=anglecoeffs[j][3]
                if (atype==btype[0]):
                    anglecoeffs.append([angletypeID,ksita,sita0,atomtype]) 
                    flag=1
                    break
            if flag==0:
                anglecoeffs.append([angletypeID,50.0,109.4710,atomtype])  
    # sorting    
    anglecoeffs.sort()

    return anglecoeffs,warning

# Get dihs coeffs: harmonic
def getDihsCoeffs(dihstypes):

    warning=""
    type_done=[]
    dihscoeffs=[]
    flag=0

    Forcefieldfile=getDreidingParamFile()
    fin=open(Forcefieldfile,'r')
    dataline=fin.readline()
    while dataline!="" and dataline!="\n" and flag==0:
        words=dataline[0:len(dataline)-2]
        if str(words).upper()=="TORSIONS":
            flag=1
            dataline=fin.readline()
            words=str.split(dataline[0:len(dataline)-1])
            while str(words[0]).upper()!="END":
                atype1=str(words[1])
                atype2=str(words[2])
                atype0=str(words[0])
                atype3=str(words[3])

                napirs=1
                if len(atype1)>2 and len(atype2)>2:
                    a1third=atype1[2]
                    a2third=atype2[2]
                    if a1third=="R":
                        a1third='2'
                    if a2third=="R":
                        a2third='2'
                    if (a1third in ['2','3']) and (a2third in ['2','3']):
                        npairs=eval(a1third)*eval(a2third)
                    else:
                        npairs=0
                        sys.exit("Error: Torsion bond has no torsion stiffness.")
                kv=0.5*1.0/npairs*eval(words[5])
                nv=eval(words[6])
                dv=(-1)*eval(words[7])            ##Lammps use a different equation for torsion,d has a opposite sign.

                for i in range(len(dihstypes)):
                    dihstypeID=dihstypes[i][0]
                    atom1type=dihstypes[i][1]
                    atom2type=dihstypes[i][2]
                    atom3type=dihstypes[i][3]
                    atom4type=dihstypes[i][4]
                    if atom2type=='O_3' or atom3type=='O_3':
                        atom1type='X'
                        atom4type='X'
                    if atom2type=='C_R' or atom3type=='C_R':
                        atom1type='X'
                        atom4type='X'
                    if atom2type=='N_R' or atom3type=='C_2':
                        atom1type='X'
                        atom4type='X'
                    if atom2type=='N_2' or atom3type=='C_2':
                        atom1type='X'
                        atom4type='X'
                    if (atype0[0:3]==atom1type[0:3] and atype1==atom2type and atype2==atom3type and atype3[0:3]==atom4type[0:3]):
                        if dihstypeID not in type_done:
                            type_done.append(dihstypeID)
                            dihscoeffs.append([dihstypeID,kv,nv,dv,atom1type,atom2type,atom3type,atom4type])
                    if (atype0[0:3]==atom4type[0:3] and atype2==atom2type and atype1==atom3type and atype3[0:3]==atom1type[0:3]):
                        if dihstypeID not in type_done:
                            type_done.append(dihstypeID)
                            dihscoeffs.append([dihstypeID,kv,nv,dv,atom1type,atom2type,atom3type,atom4type])
                dataline=fin.readline()
                words=str.split(dataline[0:len(dataline)-1])
        dataline=fin.readline()
    fin.close()
    # if not assigned
    assigned=[]
    for j in range(len(dihscoeffs)):
        dihstypeID=dihscoeffs[j][0]
        assigned.append(dihstypeID)
    for i in range(len(dihstypes)):
        dihstypeID=dihstypes[i][0]
        atom3type=dihstypes[i][3]
        atom2type=dihstypes[i][2]
        atom1type=dihstypes[i][4]
        atom4type=dihstypes[i][1]
        if dihstypeID not in assigned:
            warning=warning+" Dihstype "+str(dihstypeID)+" might be not assigned correctly."
            flag=0
            a3type=atom3type[0]
            a2type=atom2type[0]
            for j in range(len(dihscoeffs)):
                btypeID=dihscoeffs[j][0]
                kv=dihscoeffs[j][1]
                nv=dihscoeffs[j][2]
                dv=dihscoeffs[j][3]
                b3type=dihscoeffs[j][6]
                b2type=dihscoeffs[j][5]
                if (a3type==b3type[0] and a2type==b2type[0]) or (a3type==b2type[0] and a2type==b3type[0]):
                    dihscoeffs.append([dihstypeID,kv,nv,dv,atom1type,atom2type,atom3type,atom4type]) 
                    flag=1
                    break
            if flag==0:
                dihscoeffs.append([dihstypeID,0.11111,3.0,1.0,atom1type,atom2type,atom3type,atom4type])  
    # sorting    
    dihscoeffs.sort()

    return dihscoeffs,warning

def getImpsCoeffs(impstypes):
    warning=""
    impscoeffs=[]
    imps_nonzero=[]
    flag=0

    Forcefieldfile=getDreidingParamFile()
    fin=open(Forcefieldfile,'r')
    dataline=fin.readline()
    while dataline!="" and dataline!="\n" and flag==0:
        words=dataline[0:len(dataline)-2]
        if str(words).upper()=="INVERSIONS":
            flag=1
            dataline=fin.readline()
            words=str.split(dataline[0:len(dataline)-1])
            while str(words[0]).upper()!="END":
                atype=str(words[0])
                ksita=0.5*eval(words[5])
                sita0=eval(words[6])

                for i in range(len(impstypes)):
                    impstypeID=impstypes[i][0]
                    atomtype=impstypes[i][1]
                    if (atype==atomtype):
                        impscoeffs.append([impstypeID,ksita,sita0,atomtype])
                        imps_nonzero.append(impstypeID)
                dataline=fin.readline()
                words=str.split(dataline[0:len(dataline)-1])
        dataline=fin.readline()
    fin.close()
    nimps=len(impstypes)
    for i in range(nimps):
        i1=i+1
        if i1 not in imps_nonzero:
            impscoeffs.append([i1,0.0,0.0,impstypes[i1-1][1]])
    # if not assigned
    assigned=[]
    for j in range(len(impscoeffs)):
        impstypeID=impscoeffs[j][0]
        assigned.append(impstypeID)
    for i in range(len(impstypes)):
        impstypeID=impstypes[i][0]
        atomtype=impstypes[i][1]
        if impstypeID not in assigned:
            warning=warning+" Impstype "+str(impstypeID)+" might be not assigned correctly."
            flag=0
            atype=atomtype[0]
            for j in range(len(impscoeffs)):
                btypeID=impscoeffs[j][0]
                ksita=impscoeffs[j][1]
                sita0=impscoeffs[j][2]
                btype=impscoeffs[j][3]
                if (atype==btype[0]):
                    impscoeffs.append([impstypeID,ksita,sita0,atomtype]) 
                    imps_nonzero.append(impstypeID)
                    flag=1
                    break

    return impscoeffs,warning

def getPairCoeffs(atomtypes):

    warning=""
    n=len(atomtypes)
    vdws=[]
    flag=0

    Forcefieldfile=getDreidingParamFile()
    fin=open(Forcefieldfile,'r')
    dataline=fin.readline()
    while dataline!="" and dataline!="\n" and flag==0:
        words=dataline[0:len(dataline)-2]
        if str(words).upper()=="DIAGONAL_VDW":
            flag=1
            dataline=fin.readline()
            words=str.split(dataline[0:len(dataline)-1])
            while str(words[0]).upper()!="END":
                atype=str(words[0])
                R0=eval(words[2])
                D0=eval(words[3])
                csi=eval(words[4])
                for i in range(n):
                    atomtypeID=atomtypes[i][0]
                    atomtype=atomtypes[i][1]
                    if (atype==atomtype):
                        vdws.append([atomtypeID,R0,D0,csi,atomtype])
                dataline=fin.readline()
                words=str.split(dataline[0:len(dataline)-1])
        dataline=fin.readline()
    fin.close()
    
    assigned=[]
    for j in range(len(vdws)):
        atomtypeID=vdws[j][0]
        assigned.append(atomtypeID)
    for i in range(n):
        atomtypeID=atomtypes[i][0]
        atomtype=atomtypes[i][1]
        if atomtypeID not in assigned:
            warning=warning+" Atomtype "+str(atomtypeID)+" might be not assigned correctly."
    # Generate pair coeffs for input scripts
    paircoeffs=[]
    R0=[]
    D0=[]
    csi=[]
    atype=[]
    A=[]
    B=[]
    C=[]
    for i in range(n+1):
        R0.append(3.5)
        D0.append(0.1)
        csi.append(12.0)
        atype.append(atomtypes[i-1][1])
        A.append(i)
        B.append(i)
        C.append(i)

    for j in range(len(vdws)):
        atomtypeID=vdws[j][0]
        R0[atomtypeID]=vdws[j][1]
        D0[atomtypeID]=vdws[j][2]
        csi[atomtypeID]=vdws[j][3]
        atype[atomtypeID]=vdws[j][4]
    for i in range(1,n+1):
        A[i]=D0[i]*(6.0/(csi[i]-6.0))*exp(csi[i])
        B[i]=D0[i]*(csi[i]/(csi[i]-6.0))*R0[i]**6
        C[i]=csi[i]/R0[i]

    #cwd=os.path.dirname(os.path.abspath(sys.argv[0]))
    #datadir=os.path.join(cwd, '..', '..', '..','..', 'data')
    #filename=os.path.join(datadir, 'paircoeffs.txt')
    filename='X6paircoeffs.txt'
    fout=open(filename,'w')
    for i in range(1,n+1):
        for j in range(i,n+1):
            Aij=(A[i]*A[j])**0.5
            Bij=(B[i]*B[j])**0.5
            Cij=0.5*(C[i]+C[j])
            rouij=1.0/Cij
            comment=atype[i]+"  "+atype[j]
            print >>fout,"pair_coeff\t",i,"\t",j,"\t",Aij,rouij,Bij,"  #",comment
    fout.close()

    # LJ pair
    filename='LJpaircoeffs.txt'
    fout=open(filename,'w')
    for i in range(1,n+1):
        epsilon=D0[i]
        c=2**(1.0/6.0)
        sigma=R0[i]/c
        comment=atype[i]
        print >>fout,"pair_coeff\t",i,epsilon,sigma,"  #",comment
    fout.close()

    return warning
####################################################################################
