# struc2lammps

This tool generates a [data file for LAMMPS](https://lammps.sandia.gov/doc/2001/data_format.html#_cch3_930958962) simulations of molecular systems starting with an atomistic structure.

### Tool features:
- Accepts a wide range of input formats (e.g. CIF, PDB, XYZ)
- Works with periodic boundary conditions, even when covalent bonds cross boundariesThe tool uses OpenBabel to convert the initial files
- Sets up energy expression for Dreiding and PCFF force fields (including automatic bond recognition and atom typing)
- Assigns atomic changes using the Gasteiger method

*Known bug: If resonance ring crosses the boundary, atoms will be incorrectly*

## How to run

Please be aware of the **requirements** specified below, before running struc2lammps. If ready to run, then follow the instructions:

1. Specify the input arguments in [input.args](./input.args)
   - Initial structure (can be added in ./data/)
   - If box parameters are not included in the initial structure please edit the [box.txt](./box.txt) file (include cell dimensions and angles)
   - Bondscale parameter for connectivity identification (scaling sum of covalent radii of bonded atoms)
   - Force field. Currently supporting Dreiding and PCFF
2. run `python str2lmp.py`
3. Just wait! You should see a **LAMMPSDataFile.data** once if finishes.

The data file provides pair coefficents in LJ (12/6) form but you will obtain a extra file with Buckingham potential (exp/6)

## Requirements:
In order to sucesfully run **struc2lammps** and obtain a data file you must have:

- Properly installed [OpenBabel](http://openbabel.org/wiki/Category:Installation) with python bindings. Remember that struc2lammps only works with [file formats currently supported by OpenBabel](http://openbabel.org/docs/current/FileFormats/Overview.html).
  
  For more information about installing python bindings you can visit: https://openbabel.org/docs/dev/Installation/install.html#compile-bindings or https://pypi.org/project/openbabel/
 
- You also need to have installed [OVITO](https://ovito.org). This is particularly important if your structure is not orthogonal. You will be prompted to specify the location of your OVITO’s Python interpreter. Information about where to find the **ovitos** script can be found here: https://ovito.org/manual/python/introduction/running.html
