import os

def run_data4lammps(chargeMethod, infile_is_pdb, boundaries, ffname):

    data4lammps_dir = './data4lammps/bin'
    path = os.path.join(data4lammps_dir, 'main.py')
    path2 = os.path.join(data4lammps_dir, 'doAtomTyping.py')    
    new_typing_command = 'python {0} {1}'.format(path2, ffname)
    print('Generating {} atom types'.format(ffname.upper()))  
    if infile_is_pdb:
        cell_sizes = boundaries
        data4lammps_command = 'python2.7 {0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10}'.format(path, \
                    cell_sizes[0], cell_sizes[1],\
                    cell_sizes[2], cell_sizes[3],\
                    cell_sizes[4], cell_sizes[5],\
                    round(cell_sizes[6],4), round(cell_sizes[7],4), round(cell_sizes[8],4), chargeMethod)
    else:
        cell_sizes = read_cell_sizes('./bonds/bonded.lmpdat')
        if 'tilt' in cell_sizes:
            data4lammps_command = 'python2.7 {0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10}'.format(path, \
                    cell_sizes['x'][0], cell_sizes['x'][1],\
                    cell_sizes['y'][0], cell_sizes['y'][1],\
                    cell_sizes['z'][0], cell_sizes['z'][1],\
                    cell_sizes['tilt'][0], cell_sizes['tilt'][1], cell_sizes['tilt'][2], chargeMethod)
        else:        
            data4lammps_command = 'python2.7 {0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10}'.format(path, \
                    cell_sizes['x'][0], cell_sizes['x'][1],\
                    cell_sizes['y'][0], cell_sizes['y'][1],\
                    cell_sizes['z'][0], cell_sizes['z'][1],\
                    0.0, 0.0, 0.0, chargeMethod)

    return new_typing_command, data4lammps_command

def read_cell_sizes(data_file):

    cell_sizes = {}
    with open(data_file) as f:
        for line in f:

            if line.endswith('xhi\n'):
                values = line.rstrip().split()
                lo = values[0]
                hi = values[1]
                cell_sizes['x'] = [lo, hi]

            elif line.endswith('yhi\n'):
                values = line.rstrip().split()
                lo = values[0]
                hi = values[1]
                cell_sizes['y'] = [lo, hi]

            elif line.endswith('zhi\n'):
                values = line.rstrip().split()
                lo = values[0]
                hi = values[1]
                cell_sizes['z'] = [lo, hi]

            elif line.endswith('yz\n'):               
                values = line.rstrip().split()
                xy = values[0]
                xz = values[1]
                yz = values[1]         
                cell_sizes['tilt'] = [xy, xz, yz]
                
                break

    return cell_sizes

def read_cell_sizes_pdb(connected_pdb):

    cell_sizes = {}
    with open(connected_pdb) as f:
        for line in f:
            if line.startswith('CRYST1'):

                values = line.rstrip().split()[1:]
                cell_sizes['x'] = [0, values[0]]
                cell_sizes['y'] = [0, values[1]]
                cell_sizes['z'] = [0, values[2]]

                break

    return cell_sizes
