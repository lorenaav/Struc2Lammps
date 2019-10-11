import os
import string
        
###########################################################

def minCoords(atoms):

    xmin=1e10
    ymin=1e10
    zmin=1e10
    for i in range(len(atoms)):
        x=atoms[i][4]
        y=atoms[i][5]
        z=atoms[i][6]
        if x<xmin:
            xmin=x
        if y<ymin:
            ymin=y
        if z<zmin:
            zmin=z

    return xmin,ymin,zmin

def maxCoords(atoms):

    xmax=-1e10
    ymax=-1e10
    zmax=-1e10
    for i in range(len(atoms)):
        x=atoms[i][4]
        y=atoms[i][5]
        z=atoms[i][6]
        if x>xmax:
            xmax=x
        if y>ymax:
            ymax=y
        if z>zmax:
            zmax=z
    
    return xmax,ymax,zmax

# Get atomtype information
def Atomtypes(inpfile):  

    atomtypes=[]
    
    fin=open(inpfile,'r')
    dataline=fin.readline()
    while dataline!="" and dataline!="\n":
        words = string.split(dataline[0:len(dataline)-1])
        atomtypeID=eval(words[0])
        atomtype=words[1]
        if atomtype=="H__HB":
            atomtype="H___b"
        if atomtype=="H__HA":
            atomtype="H___A"
        atomtypes.append([atomtypeID,atomtype])
        dataline=fin.readline()
    atomtypes.sort()    

    return atomtypes

# Sort atom information
def AtomsInfo(inpfile):  

    fin=open(inpfile,'r')
    dataline=fin.readline()
    words = string.split(dataline)
    totalatoms=eval(words[0])
    atoms=[]
    for i in range(totalatoms):
        atoms.append([])
        
    dataline=fin.readline()
    words = string.split(dataline,sep="\n")
    while dataline!="":
        if (str(words[0]).strip()).upper()=='ATOMS':
            dataline=fin.readline()
            dataline=fin.readline()
            while dataline!="\n" and dataline!="":
                words = string.split(dataline[0:len(dataline)-1])
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
                dataline=fin.readline()
        dataline=fin.readline()
        words = string.split(dataline,sep="\n")
    fin.close()

    return atoms

# Find links to each atom
def AtomLink(atoms,bonds):   

    totalatoms=len(atoms)
    nlink=[]
    atomlinktable=[]
    for i in range(totalatoms):
        nlink.append(0)
        atomlinktable.append([0,0,0,0])
        
    for k in range(len(bonds)):
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
