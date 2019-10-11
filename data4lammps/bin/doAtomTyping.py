import os
import sys
from numpy import dot, cross

forcefield=sys.argv[1]

def getAtomtypes():

    atomtypes=[]
    fin=open("./types/atom_type.dat",'r')
    for line in fin:
        words=str.split(line[0:len(line)-1])
        atype=str(words[1]) # Dont change the upper because otherwise wont recognize several values
        #atype=str(words[1]).upper()
        if len(atype)==1:
            atype=atype+'_'
        atomtypes.append([words[0],atype])

    return atomtypes    

def getRingatoms(N):

    atominring=[]
    for i in range(N):
        atominring.append(0)
        
    fin=open("./types/atominring.dat",'r')
    for line in fin:
        words=str.split(line[0:len(line)-1])
        atomID=eval(words[0])
        ringsize=eval(words[1])
        atominring[atomID-1]=ringsize

    return atominring    

def getBondorders(NB):

    bondorder=[]
    for i in range(NB):
        bondorder.append(0)
        
    fin=open("./types/bondorder.dat",'r')
    for line in fin:
        words=str.split(line[0:len(line)-1])
        bondID=eval(words[0])
        BOvalue=eval(words[1])
        bondorder[bondID-1]=BOvalue

    return bondorder

# Find links to each atom
def AtomLink(atoms,bonds):   

    totalatoms=len(atoms)
    nlink=[]
    atomlinktable=[]
    for i in range(totalatoms):
        nlink.append(0)
        atomlinktable.append([0,0,0,0])
        
    for k in range(len(bonds)):
        #bondID=bonds[k][0]
        #bondtype=bonds[k][1]
        atom1=bonds[k][2]
        atom2=bonds[k][3]
        m1=nlink[atom1-1]
        m2=nlink[atom2-1]
        atomlinktable[atom1-1][m1]=atom2
        atomlinktable[atom2-1][m2]=atom1
        nlink[atom1-1]=nlink[atom1-1]+1
        nlink[atom2-1]=nlink[atom2-1]+1
    
    atomlinks=[]
    for k in range(len(atomlinktable)):
        atomlinks.append([])
        for m in range(4):
            atomID=atomlinktable[k][m]
            if atomID!=0:
                atomlinks[k].append(atomID)

    return atomlinks

def linkelems(atomlinks,atoms,atomtypes):

    atomlinkelems=[]
    for i in range(len(atomlinks)):
        atomlinkelems.append([])
        for j in range(len(atomlinks[i])):
            atomID=atomlinks[i][j]
            atomtypeID=atoms[atomID-1][2]
            atomtype=atomtypes[atomtypeID-1][1]
            elemtype=atomtype[0:2]
            atomlinkelems[i].append(elemtype)

    return atomlinkelems        

# Sort atom information
def getAtoms():  

    fin=open("./types/atoms.dat",'r')
    dataline=fin.readline()
    words = str.split(dataline)
    totalatoms=eval(words[0])
    atoms=[]
    for i in range(totalatoms):
        atoms.append([])
        
    dataline=fin.readline()
    dataline=fin.readline()
    dataline=fin.readline()
    for i in range(totalatoms):
        dataline=fin.readline()
        words = str.split(dataline[0:len(dataline)-1])
        natom=eval(str(words[0]))
        nmolecule=eval(str.strip(str(words[1])))
        atomtype=eval(str(words[2]))
        charge=eval(str(words[3]))
        x=eval(str(words[4]))
        y=eval(str(words[5]))
        z=eval(str(words[6]))
        if len(words)>7:
            nx=eval(str(words[7]))
            ny=eval(str(words[8]))
            nz=eval(str(words[9]))
            atoms[natom-1]=[natom,nmolecule,atomtype,charge,x,y,z,nx,ny,nz]
        else:
            atoms[natom-1]=[natom,nmolecule,atomtype,charge,x,y,z]
    fin.close()

    return atoms

def getBonds():   

    fin=open("./types/bonds.dat",'r')
    bonds=[]
        
    dataline=fin.readline()
    dataline=fin.readline()
    
    bondID=0
    dataline=fin.readline()
    words = str.split(dataline[0:len(dataline)-1])
    while words!=[]:
        bondtype=eval(str(words[1]))
        atom1=eval(str(words[2]))
        atom2=eval(str(words[3]))
        
        bondID=bondID+1
        bonds.append([bondID,bondtype,atom1,atom2])
            
        dataline=fin.readline()
        words = str.split(dataline[0:len(dataline)-1])

    dataline=fin.readline()    
    words = str.split("\n")
    fin.close()

    return bonds

def getDREIDINGatomtypes(atoms,atomlinks,atomelems,atomlinkelems,atominring):

    atomtype=[]
    for i in range(len(atomelems)):
        atomtype.append('')
    for i in range(len(atomelems)):
        atomelem=atomelems[i]
        if atomelem=="H_":
            if ("O_" in atomlinkelems[i] or "N_" in atomlinkelems[i] or "S_" in atomlinkelems[i]): 
                atomtype[i]='H__HB'
            else:
                atomtype[i]='H_'                
        elif atomelem=="C_":
            if len(atomlinks[i])==4: 
                atomtype[i]='C_3'
            elif len(atomlinks[i])==3: 
                if atominring[i]<0:
                    atomtype[i]='C_R'                    
                else:
                    atomtype[i]='C_2'
            elif len(atomlinks[i])==2: 
                atomtype[i]='C_1'
            else:    
                atomtype[i]='C_3'                
        elif atomelem=="N_":
            if len(atomlinks[i])<=2:  
                atomtype[i]='N_1'
            elif len(atomlinks[i])==3:
                atomID0=i+1
                atomID1=atomlinks[i][0]
                atomID2=atomlinks[i][1]
                atomID3=atomlinks[i][2]
                v1=[(atoms[atomID1][4]-atoms[atomID0][4]),(atoms[atomID1][5]-atoms[atomID0][5]),(atoms[atomID1][6]-atoms[atomID0][6])]
                v2=[(atoms[atomID2][4]-atoms[atomID0][4]),(atoms[atomID2][5]-atoms[atomID0][5]),(atoms[atomID2][6]-atoms[atomID0][6])]
                v3=[(atoms[atomID3][4]-atoms[atomID0][4]),(atoms[atomID3][5]-atoms[atomID0][5]),(atoms[atomID3][6]-atoms[atomID0][6])]
                vx=cross(v2,v3)
                vp=dot(v1,vx)
                if abs(vp)<1.0e-2:  
                    atomtype[i]='N_2'
                elif atominring[i]<0:  
                    atomtype[i]='N_R'
                else:  
                    atomtype[i]='N_3'
            elif len(atomlinks[i])==4:
                atomtype[i]='N_3'
            else:
                atomtype[i]='N_3'                
        elif atomelem=="O_":
            if len(atomlinks[i])==1: 
                atomtype[i]='O_2'                    
            elif len(atomlinks[i])==2: 
                if atominring[i]<0: 
                    atomtype[i]='O_R'   
                else:    
                    atomtype[i]='O_3'                   
            elif len(atomlinks[i])==3:
                atomtype[i]='O_3' 
            else:    
                atomtype[i]='O_3' 
        elif atomelem=="S_":
            if len(atomlinks[i])==4:
                atomtype[i]="S_3"                                        
            else:
                atomtype[i]="S_"  
        elif atomelem=="P_":
            if len(atomlinks[i])==4:
                atomtype[i]="P_3"  
            else:
                atomtype[i]="P_"  
        elif atomelem=="F_":
            atomtype[i]="F_"                                        
        elif atomelem=="CL":
            atomtype[i]="Cl"                                        
        elif atomelem=="B_":
            if len(atomlinks[i])==4:
                atomtype[i]="B_3"  
            elif len(atomlinks[i])==3:
                atomtype[i]="B_2"  
            else:    
                atomtype[i]="B_3"
        else:
            atomtype[i] = atomelem
                
    return atomtype

def getPCFFatomtypes(atomlinks,atomelems,atomlinkelems,atominring,bondorder):

    atomPCFFtype=[]
    for i in range(len(atomelems)):
        atomPCFFtype.append('')
    for i in range(len(atomelems)):
        atomelem=atomelems[i]
        if atomelem=="H_":
            if ("C_" in atomlinkelems[i] or "SI" in atomlinkelems[i] or "H_" in atomlinkelems[i]): 
                atomPCFFtype[i]='h'
            elif ("O_" in atomlinkelems[i] or "N_" in atomlinkelems[i]): 
                atomPCFFtype[i]='h*'
            elif  ("AL" in atomlinkelems[atomlinks[i][0]]):
                atomPCFFtype[i]='hoa'
            elif  ("SI" in atomlinkelems[atomlinks[i][0]]):
                atomPCFFtype[i]='hos'
            else:
                atomPCFFtype[i]='h'                
        elif atomelem=="C_":
            if len(atomlinks[i])==4: 
                atomPCFFtype[i]='c'
            elif len(atomlinks[i])==3: 
                if atomlinkelems[i].count('N_')==3:
                    atomPCFFtype[i]='cr'                    
                elif atomlinkelems[i].count('O_')==3:
                    atomPCFFtype[i]='cz'                    
                elif atomlinkelems[i].count('O_')==2 and atomlinkelems[i].count('N_')==1:
                    atomPCFFtype[i]='c_2'                    
                elif atomlinkelems[i].count('O_')==1 and atomlinkelems[i].count('N_')==2:
                    atomPCFFtype[i]='c_2'                    
                elif atomlinkelems[i].count('O_')==2 and atomlinkelems[i].count('C_')==1:
                    atomPCFFtype[i]='c_1'                    
                elif atomlinkelems[i].count('O_')==1 and atomlinkelems[i].count('N_')==1 and atomlinkelems[i].count('C_')==1:
                    atomPCFFtype[i]='c_1'                    
                elif atomlinkelems[i].count('O_')==1 and atomlinkelems[i].count('C_')==2:
                    atomPCFFtype[i]='c_0'                    
                elif atomlinkelems[i].count('O_')==1 and atomlinkelems[i].count('H_')==1 and atomlinkelems[i].count('C_')==1:
                    atomPCFFtype[i]='c_0'                    
                else:
                    atomPCFFtype[i]='cp'
            elif len(atomlinks[i])==2: 
                atomPCFFtype[i]='ct'                
            if atomPCFFtype[i]=='':
                atomPCFFtype[i]='c'                                            
        elif atomelem=="N_":
            if atominring[i]==5 or atominring[i]==6:  # in a 5- or 6-ring
                atomPCFFtype[i]='np'
            elif atominring[i]==6 and atomlinkelems[i].count('H_')==1:  # in a 6-ring
                atomPCFFtype[i]='nh'
            if len(atomlinks[i])==1:  # sp triple
                atomPCFFtype[i]='nt'
            if len(atomlinks[i])==4: 
                atomPCFFtype[i]='n+'
            if len(atomlinks[i])==2 and atominring[i]==0:  # one double bond
                atomPCFFtype[i]='n=2'
            elif len(atomlinks[i])==2:  
                atomPCFFtype[i]='nz'
            if len(atomlinks[i])==3: 
                if atomlinkelems[i].count('H_')==2:
                    if atomlinkelems[i][0]!="H_":
                        CatomID=atomlinks[i][0]
                    elif atomlinkelems[i][1]!="H_":
                        CatomID=atomlinks[i][1]
                    elif atomlinkelems[i][2]!="H_":
                        CatomID=atomlinks[i][2]
                    if atominring[CatomID-1]>0:
                        atomPCFFtype[i]='nb'
                elif atomlinkelems[i].count('C_')>=2:
                    for k in range(3):
                        CatomID=atomlinks[i][k]
                        if atominring[CatomID-1]>0:
                            atomPCFFtype[i]='nn'                        
                if atomlinkelems[i].count('H_')==2 and atomlinkelems[i].count('C_')==1:
                    if atomlinkelems[i][0]=="C_":
                        CatomID=atomlinks[i][0]
                    elif atomlinkelems[i][1]=="C_":   
                        CatomID=atomlinks[i][1]                        
                    elif atomlinkelems[i][2]=="C_":   
                        CatomID=atomlinks[i][2]                        
                    if len(atomlinks[CatomID-1])==3 and atomlinkelems[CatomID-1].count('N_')==3:
                        atomPCFFtype[i]='nr'
                    else:
                        atomPCFFtype[i]='na'                        
            if atomPCFFtype[i]=='':
                atomPCFFtype[i]='n'                                            
        elif atomelem=="O_":
            if len(atomlinks[i])==1: 
                if atomlinkelems[i].count('O_')==1 or atomlinkelems[i].count('N_')==1:
                    atomPCFFtype[i]='o='                    
                elif atomlinkelems[i].count('S_')==1 or atomlinkelems[i].count('P_')==1:
                    atomPCFFtype[i]='o='                    
                elif atomlinkelems[i].count('C_')==1:
                    CatomID=atomlinks[i][0]
                    if len(atomlinks[CatomID-1])==3:
                        atomPCFFtype[i]='o_1'                    
                    else:    
                        atomPCFFtype[i]='o='                    
            if len(atomlinks[i])==2: 
                if atomlinkelems[i].count('H_')==1:
                    atomPCFFtype[i]='oh'                    
                if atomlinkelems[i].count('H_')==2:
                    atomPCFFtype[i]='o*'                    
                if atomlinkelems[i].count('C_')==2:
                    atomPCFFtype[i]='op'                    
                elif atomlinkelems[i].count('SI')==2:
                    atomPCFFtype[i]='oss'                    
                elif atomlinkelems[i].count('SI')==1 and atomlinkelems[i].count('H_')==1:
                    atomPCFFtype[i]='osh'                    
                elif atomlinkelems[i].count('AL')==1 and atomlinkelems[i].count('H_')==1:
                    atomPCFFtype[i]='oah'                    
                elif atomlinkelems[i].count('AL')==1 and atomlinkelems[i].count('SI')==1:
                    atomPCFFtype[i]='oas'                    
                elif atomlinkelems[i].count('C_')==1:
                    if atomlinkelems[i][0]=="C_":
                        CatomID=atomlinks[i][0]
                    elif atomlinkelems[i][1]=="C_":   
                        CatomID=atomlinks[i][1]                        
                    if len(atomlinks[CatomID-1])==3 and atomlinkelems[CatomID-1].count('O_')==2:
                        atomPCFFtype[i]='o_2'                    
                elif atomlinkelems[i].count('C_')==1:
                    if atomlinkelems[i][0]=="C_":
                        CatomID=atomlinks[i][0]
                    elif atomlinkelems[i][1]=="C_":   
                        CatomID=atomlinks[i][1]                        
                    if len(atomlinks[CatomID-1])==3 and atomlinkelems[CatomID-1].count('O_')==3:
                        atomPCFFtype[i]='oz'                    
                else:
                    atomPCFFtype[i]='o'                     
            if atomPCFFtype[i]=='':
                atomPCFFtype[i]='o'                                            
        elif atomelem=="S_":
            if len(atomlinks[i])==1 and atomlinkelems[i].count('C_')==1:
                    atomPCFFtype[i]="s'"                                        
            if len(atomlinks[i])==2: 
                if atomlinkelems[i].count('C_')==2:
                    atomPCFFtype[i]='sc'                    
                elif atomlinkelems[i].count('C_')==2: # another rule
                    atomPCFFtype[i]='sp'                    
                elif atomlinkelems[i].count('H_')==1:
                    atomPCFFtype[i]='sh'                    
                else:     
                    atomPCFFtype[i]='s'                    
            if len(atomlinks[i])==4: 
                if atomlinkelems[i].count('O_')==1:
                    atomPCFFtype[i]="sf"     
            if atomPCFFtype[i]=='':
                atomPCFFtype[i]='s'                                            
        elif atomelem=="P_":
            if atomlinkelems[i].count('N_')==4: 
                atomPCFFtype[i]='p='
            else: 
                atomPCFFtype[i]='p'
            if atomPCFFtype[i]=='':
                atomPCFFtype[i]='p'                                            
        #elif atomelem=="F_":
        #    atomPCFFtype[i]="f"                                        
        #elif atomelem=="CL":
        #    atomPCFFtype[i]="cl"                                        
        else:
            atomPCFFtype[i]=remove_underscore(atomelem)

    return atomPCFFtype 

def collectAtomtypes(atomtypes):

    sumatomtypes=[]
    for i in range(len(atomtypes)):
        atomtype=atomtypes[i]
        if atomtype not in sumatomtypes:
            sumatomtypes.append(atomtype)    

    return sumatomtypes        

def collectBondtypes(bondtypes):

    sumbondtypes=[]
    for i in range(len(bondtypes)):
        i1=bondtypes[i][0]
        i2=bondtypes[i][1]
        
        bondtype1=[i1,i2]
        bondtype2=[i2,i1]
        if bondtype1 not in sumbondtypes and bondtype2 not in sumbondtypes:
            sumbondtypes.append(bondtype1)    

    return sumbondtypes        

def elems(atoms,atomtypes):

    atomelems=[]
    for i in range(len(atoms)):
        atomtypeID=atoms[i][2]
        atomtype=atomtypes[atomtypeID-1][1]
        elem=atomtype[0:2]
        atomelems.append(elem)

    return atomelems

def assignAtomtypeID(newatomtypes,sumatomtypes):

    atomwithtypeID=[]
    for i in range(len(newatomtypes)):
        atomwithtypeID.append(0)
    for k in range(len(sumatomtypes)):
        typeID=k+1
        stype=sumatomtypes[k]
        for i in range(len(newatomtypes)):
            atype=newatomtypes[i]
            if atype==stype:
                atomwithtypeID[i]=typeID

    return atomwithtypeID

def assignBondtypeID(newbondtypes,atoms,bonds):

    for k in range(len(newbondtypes)):
        typeID=k+1
        k1=newbondtypes[k][0]
        k2=newbondtypes[k][1]
        stype=[k1,k2]
        for i in range(len(bonds)):
            i1=bonds[i][2]
            i2=bonds[i][3]
            a1type=atoms[i1-1][2]
            a2type=atoms[i2-1][2]
            btype1=[a1type,a2type]
            btype2=[a2type,a1type]
            if btype1==stype or btype2==stype:
                bonds[i][1]=typeID

    return bonds

def atomTyping(forcefield,atoms,atomlinks,atomelems,atomlinkelems,atominring,bondorder):

    if forcefield=="PCFF":
        newatomtypes=getPCFFatomtypes(atomlinks,atomelems,atomlinkelems,atominring,bondorder)    
    else:
        newatomtypes=getDREIDINGatomtypes(atoms,atomlinks,atomelems,atomlinkelems,atominring)   
        
    # summarize    
    sumatomtypes=collectAtomtypes(newatomtypes)    
    atomwithnewtypeID=assignAtomtypeID(newatomtypes,sumatomtypes)

    return sumatomtypes,atomwithnewtypeID

def updateAtoms(atoms,atomwithnewtypeID):

    for i in range(len(atoms)):
        atoms[i][2]=atomwithnewtypeID[i]

    return atoms

def bondTyping(atomwithnewtypeID,atoms,bonds):

    bondtypes=[]
    for i in range(len(bonds)):
        bondtypes.append([])
    for i in range(len(bonds)):
        bondID=bonds[i][0]
        i1=bonds[i][2]
        i2=bonds[i][3]
        a1type=atoms[i1-1][2]
        a2type=atoms[i2-1][2]
        bondtypes[bondID-1]=[a1type,a2type]
    # collect types    
    newbondtypes=collectBondtypes(bondtypes)
    bonds=assignBondtypeID(newbondtypes,atoms,bonds)

    return newbondtypes,bonds

def outputTyping(forcefield,atomtypes,bondtypes,atoms,bonds):
    
    with open("forcefield.name",'w') as f:
        f.write(forcefield+'\n')
                
    with open("./types/newatom_type.dat",'w') as fout:
        
        assert len(atomtypes) > 0
        for i in range(len(atomtypes)):
            ntype=i+1
            atype=atomtypes[i]
            fout.write('{0} {1}\n'.format(ntype,atype))
    
    with open("./types/newbond_type.dat",'w') as fout:
        for i in range(len(bondtypes)):
            ntype=i+1
            a1typeID=bondtypes[i][0]
            a2typeID=bondtypes[i][1]
            a1type=atomtypes[a1typeID-1]
            a2type=atomtypes[a2typeID-1]        
            fout.write('{0} {1} {2}\n'.format(ntype,a1type,a2type))
    
    outfile="./types/newatoms.dat"
    fout=open(outfile,'w')
    fout.write("{0} Atoms\n".format(len(atoms)))
    fout.write("\n")
    fout.write("Atoms\n")
    fout.write("\n")
    for i in range(len(atoms)):
        outline=""
        for j in range(len(atoms[i])):
            outline=outline+str(atoms[i][j])+" "
        fout.write(outline+'\n')
    fout.close()        

    outfile="./types/newbonds.dat"
    fout=open(outfile,'w')
    fout.write("Bonds\n")
    fout.write("\n")
    for i in range(len(bonds)):
        outline=""
        for j in range(len(bonds[i])):
            outline=outline+str(bonds[i][j])+" "
        fout.write(outline+'\n')
    fout.close()
    
def remove_underscore(atomtype):

    if atomtype.endswith('_'):
        return atomtype[:-1].lower()

    return atomtype.lower()

def main(forcefield):

    atomtypes=getAtomtypes()
    atoms=getAtoms()
    bonds=getBonds()
    N=len(atoms)
    NB=len(bonds)
    atominring=getRingatoms(N)
    bondorder=getBondorders(NB)
    
    atomlinks=AtomLink(atoms,bonds)
    atomelems=elems(atoms,atomtypes)
    atomlinkelems=linkelems(atomlinks,atoms,atomtypes)
    
    newatomtypes,atomwithnewtypeID=atomTyping(forcefield,atoms,atomlinks,atomelems,atomlinkelems,atominring,bondorder)
    atoms=updateAtoms(atoms,atomwithnewtypeID)
    newbondtypes,bonds=bondTyping(atomwithnewtypeID,atoms,bonds)

    outputTyping(forcefield,newatomtypes,newbondtypes,atoms,bonds)

main(forcefield)
