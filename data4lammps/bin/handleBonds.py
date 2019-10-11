import os
import string

def outputtypes(title,typelist):

    filename=str(title)+".index"
    fout=open(filename,'w')
    for i in range(len(typelist)):
        outline=" "
        for j in range(len(typelist[i])):
            outline=outline+str(typelist[i][j])+" "
        print >>fout,i+1,outline
    fout.close()        

def createBonds(atoms,bondlength):

    totalatoms=len(atoms)
    fout=open('atoms.data','w')
    print >>fout,"Totalatoms: ",totalatoms
    print >>fout,"Bondlength: ",bondlength
    for i in range(totalatoms):
        k=atoms[i][0]
        x=atoms[i][4]
        y=atoms[i][5]
        z=atoms[i][6]
        print >>fout, '%8d %15.8e %15.8e %15.8e' % (k,x,y,z)
    fout.close()
    os.system("./findBonds")  #read 'atoms.data' and generate 'bonds.data'
    
    return 0 #nbonds

def createAngles(atomlinks):

    totalatoms=len(atomlinks)
    nangles=0
    angles=[]
    for i in range(totalatoms):
        nlinked=len(atomlinks[i])
        if nlinked>2:
            for j in range(nlinked-1):
                kl=atomlinks[i][j]
                ko=i+1
                for j1 in range(j+1,nlinked):
                    kr=atomlinks[i][j1]
                    angletype=1
                    nangles=nangles+1
                    angles.append([nangles,angletype,kl,ko,kr])
        elif nlinked==2:
            nangles=nangles+1
            kl=atomlinks[i][0]
            ko=i+1
            kr=atomlinks[i][1]
            angletype=1
            angles.append([nangles,angletype,kl,ko,kr])
            
    return angles

def createDihedrals(atomlinks,bonds):

    totalbonds=len(bonds)
    ndihs=0
    dihs=[]
    for i in range(totalbonds):
        bondID=bonds[i][0]
        bondtype=bonds[i][1]
        atom1=bonds[i][2]
        atom2=bonds[i][3]
        n1=len(atomlinks[atom1-1])
        n2=len(atomlinks[atom2-1])
        for j1 in range(n1):
            k1=atomlinks[atom1-1][j1]
            if k1!=atom2:
                k2=atom1
                k3=atom2
                for j2 in range(n2):
                    k4=atomlinks[atom2-1][j2]
                    if k4!=atom1:
                        ndihs=ndihs+1
                        dihtype=1
                        dihs.append([ndihs,dihtype,k1,k2,k3,k4])
            
    return dihs

def createImpropers(atomlinks):

    totalatoms=len(atomlinks)
    nimps=0
    imps=[]
    for i in range(totalatoms):
        nlinked=len(atomlinks[i])
        if nlinked==3:
            nimps=nimps+1
            k1=i+1
            k2=atomlinks[i][0]
            k3=atomlinks[i][1]
            k4=atomlinks[i][2]
            imptype=1
            imps.append([nimps,imptype,k1,k2,k3,k4])
        if nlinked>3:
            for j in range(nlinked):
                nimps=nimps+1
                k1=i+1
                k2=atomlinks[i][j]
                
                j1=j+1
                if j1>=nlinked:
                    j1=0
                j2=j1+1
                if j2>=nlinked:
                    j2=0                
                k3=atomlinks[i][j1]
                k4=atomlinks[i][j2]
                imptype=1
                imps.append([nimps,imptype,k1,k2,k3,k4])

    return imps

def getBondtypes(inpfile):   

    fin=open(inpfile,'r')
    bondtypes=[]
        
    dataline=fin.readline()
    words = string.split(dataline,sep="\n")
    while dataline!="\n" and dataline!="":
        words = string.split(dataline[0:len(dataline)-1])
        bondtypeID=eval(str(words[0]))
        atom1type=str(words[1])
        atom2type=str(words[2])        
        bondtypes.append([bondtypeID,atom1type,atom2type])            
        dataline=fin.readline()
    fin.close()

    return bondtypes

def getBonds(inpfile,fout,getorout):   

    fin=open(inpfile,'r')
    bonds=[]
        
    dataline=fin.readline()
    words = string.split(dataline,sep="\n")
    while str(words[0]).upper()!="ANGLES" and dataline!="":
        if str(words[0]).upper()=='BONDS':
            dataline=fin.readline()
            dataline=fin.readline()
            bondID=0
            words = string.split(dataline[0:len(dataline)-1])
            while words!=[]:
                bondtype=eval(str(words[1]))
                atom1=eval(str(words[2]))
                atom2=eval(str(words[3]))
                
                bondID=bondID+1
                if getorout==1:
                    bonds.append([bondID,bondtype,atom1,atom2])
                else:
                    print >>fout,'%8d %8d %8d %8d' % (bondID,bondtype,atom1,atom2)
                    
                dataline=fin.readline()
                words = string.split(dataline[0:len(dataline)-1])
    
        dataline=fin.readline()    
        words = string.split(dataline,sep="\n")
    fin.close()
    if getorout==1:

        return bonds

    else:

        return 0

# Read angles from datafile
def getAngles(inpfile,fout,getorout):   

    fin=open(inpfile,'r')
    angles=[]
        
    dataline=fin.readline()
    words = string.split(dataline,sep="\n")
    while str(words[0]).upper()!="DIHEDRALS" and dataline!="":
        if str(words[0]).upper()=='ANGLES':
            dataline=fin.readline()
            dataline=fin.readline()
            angleID=0
            words = string.split(dataline[0:len(dataline)-1])
            while words!=[]:
                angletype=eval(str(words[1]))
                atom1=eval(str(words[2]))
                atom2=eval(str(words[3]))
                atom3=eval(str(words[4]))
                
                angleID=angleID+1
                if getorout==1:  ##get
                    angles.append([angleID,angletype,atom1,atom2,atom3])
                else:
                    print >>fout,'%8d %8d %8d %8d %8d' % (angleID,angletype,atom1,atom2,atom3)
                dataline=fin.readline()
                words = string.split(dataline[0:len(dataline)-1])
    
        dataline=fin.readline()    
        words = string.split(dataline,sep="\n")
    fin.close()
    if getorout==1:

        return angles

    else:

        return 0

def getAngletypes(angles,atoms,atomtypes):

    angletype=[]
    angletypes=[]
    angletypeID=0
    
    nangle=len(angles)
    for i in range(nangle):
        catomID=angles[i][3]
        atomtypeID=atoms[catomID-1][2]
        if atomtypeID not in angletype:
            angletypeID=angletypeID+1
            angletype.append(atomtypeID)
            angletypes.append([angletypeID,atomtypeID])
    outputtypes("angletypes",angletypes)

    for i in range(len(angletypes)):
        atomtypeID=angletypes[i][1]
        for j in range(len(atomtypes)):
            atomtypeID0=atomtypes[j][0]
            if atomtypeID==atomtypeID0:
                atomtype=atomtypes[j][1]
        angletypes[i][1]=atomtype  
    # update angletype in angles    
    for i in range(nangle):
        catomID=angles[i][3]
        atomtypeID=atoms[catomID-1][2]
        atomtype=atomtypes[atomtypeID-1][1]
        for j in range(len(angletypes)):
            atomtype0=angletypes[j][1]
            if atomtype==atomtype0:
                angletypeID=angletypes[j][0]
                
        angles[i][1]=angletypeID        
        
    return angletypes,angles

def getDihs(inpfile,fout,getorout):   

    fin=open(inpfile,'r')
    dihs=[]
        
    dataline=fin.readline()
    words = string.split(dataline,sep="\n")
    while str(words[0]).upper()!="IMPROPERS" and dataline!="":
        if str(words[0]).upper()=='DIHEDRALS':
            dataline=fin.readline()
            dataline=fin.readline()
            dihID=0
            words = string.split(dataline[0:len(dataline)-1])
            while words!=[]:
                dihtype=eval(str(words[1]))
                atom1=eval(str(words[2]))
                atom2=eval(str(words[3]))
                atom3=eval(str(words[4]))
                atom4=eval(str(words[5]))
                
                dihID=dihID+1
                if getorout==1:
                    dihs.append([dihID,dihtype,atom1,atom2,atom3,atom4])
                else:
                    print >>fout,'%8d %8d %8d %8d %8d %8d' % (dihID,dihtype,atom1,atom2,atom3,atom4)                    
                dataline=fin.readline()
                words = string.split(dataline[0:len(dataline)-1])
    
        dataline=fin.readline()    
        words = string.split(dataline,sep="\n")
    fin.close()
    if getorout==1:

        return dihs

    else:

        return 0

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
        if (atom1type[2:3]!='2' and atom1type[2:3]!='R') or (atom2type[2:3] not in ['2','R']):
            atom1type='X'
        if (atom4type[2:3]!='2' and atom4type[2:3]!='R') or (atom3type[2:3] not in ['2','R']):
            atom4type='X'
            
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
        if (atom3type[2:3]!='2' and atom3type[2:3]!='R') or (atom1type[2:3] not in ['2','R']):
            atom3type='X'
        if (atom4type[2:3]!='2' and atom4type[2:3]!='R') or (atom2type[2:3] not in ['2','R']):
            atom4type='X'
        
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

def getImps(inpfile,fout,getorout):   

    fin=open(inpfile,'r')
    imps=[]
        
    dataline=fin.readline()
    words = string.split(dataline,sep="\n")
    while str(words[0]).upper()!="PAIR COEFFS" and dataline!="":
        if str(words[0]).upper()=='IMPROPERS':
            dataline=fin.readline()
            dataline=fin.readline()
            impID=0
            words = string.split(dataline[0:len(dataline)-1])
            while words!=[]:
                imptype=eval(str(words[1]))
                atom1=eval(str(words[2]))
                atom2=eval(str(words[3]))
                atom3=eval(str(words[4]))
                atom4=eval(str(words[5]))
                
                impID=impID+1
                if getorout==1:
                    imps.append([impID,imptype,atom1,atom2,atom3,atom4])
                else:
                    print >>fout,'%8d %8d %8d %8d %8d %8d' % (impID,imptype,atom1,atom2,atom3,atom4)
                    
                dataline=fin.readline()
                words = string.split(dataline[0:len(dataline)-1])
    
        dataline=fin.readline()    
        words = string.split(dataline,sep="\n")
    fin.close()
    if getorout:

        return imps

    else:

        return 0

def getImpstypes(imps,atoms,atomtypes):

    impstype=[]
    impstypes=[]
    impstypeID=0
    
    nimps=len(imps)
    for i in range(nimps):
        catomID=imps[i][2]
        atomtypeID=atoms[catomID-1][2]
        if atomtypeID not in impstype:
            impstypeID=impstypeID+1
            impstype.append(atomtypeID)
            impstypes.append([impstypeID,atomtypeID])

    for i in range(len(impstypes)):
        atomtypeID=impstypes[i][1]
        for j in range(len(atomtypes)):
            atomtypeID0=atomtypes[j][0]
            if atomtypeID==atomtypeID0:
                atomtype=atomtypes[j][1]
        impstypes[i][1]=atomtype  
    # update impstype in impss    
    for i in range(nimps):
        catomID=imps[i][2]
        atomtypeID=atoms[catomID-1][2]
        atomtype=atomtypes[atomtypeID-1][1]
        for j in range(len(impstypes)):
            atomtype0=impstypes[j][1]
            if atomtype==atomtype0:
                impstypeID=impstypes[j][0]
        imps[i][1]=impstypeID        
    outputtypes("impropertypes",impstypes)        

    return impstypes,imps
