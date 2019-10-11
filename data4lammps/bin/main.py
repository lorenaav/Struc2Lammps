import data4Lammps
import sys

if len(sys.argv) != 11:
    sys.stderr.write("Usage: python main.py xlo xhi ylo yhi zlo zhi xy xz yz chargeMethod\n")
    sys.exit()
forcefield=data4Lammps.getForcefield()
structureName=data4Lammps.getFileName()

if str(forcefield).upper()=="REAXFF":
    datafile=data4Lammps.createReaxDatafile(forcefield, structureName, sys.argv[1] , sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5] , sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10])

else:
    datafile=data4Lammps.createDatafile(forcefield, structureName, sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5] , sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10])
