import os
import sys
import string
import handleAtoms
import handleBonds

# Get parameter file name, with installed path
def getGasteigerParamFile():

    cwd=os.path.dirname(os.path.abspath(sys.argv[0]))
    datadir=os.path.join(cwd, '..', 'data')

    return os.path.join(datadir, 'Gasteiger_parameters.txt')
    #return 'Gasteiger_parameters.txt'

def getGasteiger_parameters(forcefield):

    Gparas=[]
    Gasteigerfile=getGasteigerParamFile()
    fin=open(Gasteigerfile,'r')
    dataline=fin.readline()
    while dataline!="":
        words=string.split(dataline[0:len(dataline)-1])
        if str(words[0]).upper()==str(forcefield).upper():
            dataline=fin.readline()
            dataline=fin.readline()
            words=string.split(dataline[0:len(dataline)-1])
            while str(words[0]).upper()!="END":
                atomtype=words[0]
                ai=eval(words[1])
                bi=eval(words[2])
                ci=eval(words[3])
                Gparas.append([atomtype,ai,bi,ci])
                dataline=fin.readline()
                words=string.split(dataline[0:len(dataline)-1])
        dataline=fin.readline()
    fin.close()

    return Gparas

def getAtomsBonds():

    inpfile='atom_type.dat'
    atomtypes=handleAtoms.Atomtypes(inpfile)
    #print "Atomtypes total=",len(atomtypes)
    inpfile='atoms.dat'
    atoms=handleAtoms.AtomsInfo(inpfile)
    #print "Atoms total=",len(atoms)       
    inpfile='bonds.dat'
    bonds=handleBonds.getBonds(inpfile,0,1)
    #print "Bonds total=",len(bonds)

    return atomtypes,atoms,bonds

def getGasteigerCharge(Gparas,atomtypes,atoms,bonds):

    natoms=len(atoms)
    nbonds=len(bonds)
    
    Q=[0]
    qn=[0]
    x=[0]
    xp=[0]
    
    # Initialization
    for i in range(natoms):
        atypeID=atoms[i][2]
        atype=atomtypes[atypeID-1][1]
        atype=atype.upper()
        ai=0
        bi=0
        ci=0
        
        found=0
        for j in range(len(Gparas)):
            Gatype=Gparas[j][0]
            Gatype=Gatype.upper()
            if atype==Gatype:
                ai=Gparas[j][1]
                bi=Gparas[j][2]
                ci=Gparas[j][3]
                found=1
                break
        if found==0:    
            atype2=atype[0]
            for j in range(len(Gparas)):
                Gatype=Gparas[j][0]
                Gatype=Gatype.upper()
                Gatype2=Gatype[0]
                if atype2==Gatype2:
                    ai=Gparas[j][1]
                    bi=Gparas[j][2]
                    ci=Gparas[j][3]
                    found=1
                    break
            
        Q.append(0)
        qn.append(0)
        x.append(ai)
        xp.append([ai,bi,ci])
        
    # Loop over n
    for n in range(1,100):
        for k in range(nbonds):
            i=bonds[k][2]
            j=bonds[k][3]
            itypeID=atoms[i-1][2]
            itype=atomtypes[itypeID-1][1]
            jtypeID=atoms[j-1][2]
            jtype=atomtypes[jtypeID-1][1]

            d1=xp[i][0]+xp[i][1]+xp[i][2]
            d2=xp[j][0]+xp[j][1]+xp[j][2]
            d=min(d1,d2)
            if itype=='H_' or jtype=='H_' or itype=='h' or jtype=='h*':
                d=max(d1,d2)
            if itype=='H___b' or jtype=='H___b':
                d=max(d1,d2)
            if itype=='H___A' or jtype=='H___A':
                d=max(d1,d2)
            if d==0:
                d=1e-6
            qn[i]=qn[i]+((x[j]-x[i])/d)*(0.5)**n
            qn[j]=qn[j]+((x[i]-x[j])/d)*(0.5)**n

        nconverg=0
        for i in range(1,natoms+1):
            if abs(qn[i])<1e-8:
                nconverg=nconverg+1
            Q[i]=Q[i]+qn[i]
            qn[i]=0.0
            x[i]=xp[i][0]+xp[i][1]*Q[i]+xp[i][2]*Q[i]*Q[i]
            
        if nconverg==natoms:
            Goutput(Q)
            break
    # Check total charge
    s=0.0
    for i in range (len(Q)):
        s=s+Q[i]
    if abs(s)>1e-10:
        print "Total charge of the system=",s
        r=s/len(Q)
        for i in range(len(Q)):
            Q[i]=Q[i]-r
        print "Total charge of the system neutralized,=0.0!"
    else:
        print "Total charge=",s

    return Q

def Goutput(Q):

    fout=open("Gasteiger_Q.txt",'w')
    print >>fout,"atom_i \t Qi"
    for i in range(1,len(Q)):
        print >>fout,'{0:18d} {1:18.8f}'.format(i,Q[i])
    fout.close()    

def main():

    Gparas=getGasteiger_parameters()
    atomtypes,atoms,bonds=getAtomsBonds()
    Q=getGasteigerCharge(Gparas,atomtypes,atoms,bonds)           
    Goutput(Q)            
#main()
       
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
