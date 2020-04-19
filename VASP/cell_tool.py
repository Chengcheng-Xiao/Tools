
#!/usr/bin/env python
"""
A script to make supercell and find atomic distance

Depends on ase and spglib
"""

import argparse
import spglib
from ase.io import read, write
import ase.io.vasp
from ase.build.supercells import make_supercell
from ase import Atoms
from ase.visualize import view
import numpy as np
import os
import sys
import datetime
import time

# is digit? function
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        pass
    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass
    return False

# Command line praser
#----------------------------
parser = argparse.ArgumentParser(description='A script to make supercell(transformed supercell) and measure atomic distances.')
parser.add_argument('-i','--insert', action="store_true", default=False, dest='over',
                    help='Overwrite POSCAR?. (Default=False)')
parser.add_argument('-vi','--visualize', action="store_true", default=False, dest="visualize",
                    help='Use ASE-GUI to visualize structure (Default=False)')
parser.add_argument('-v', '--verbose', action="store_true", default=False, dest="prt",
                    help='print out info? (Default=False)')
parser.add_argument('-c','--POSCAR', action="store", default=str('POSCAR'), dest="POSCAR",
                    help='Input file name (Default=POSCAR)')
parser.add_argument('-d', action="store_true", default=False, dest="get_super",
                    help='Make supercell? (Default=False)')
parser.add_argument('--dim', action="store", default="1 0 0 0 1 0 0 0 1", dest="T_mat",
                    help='Transformation matrix ("1,0,0,0,1,0,0,0,1" or "1 1 1")')
parser.add_argument('-di', action="store_true", default=False, dest="getdist",
                    help='Get distance? (Default=False)')
parser.add_argument('--dis', action="store", default="1 2", dest="atom_list",
                    help='Distance between atoms a and b, starting from 1 (Default="1 2")')


prm = parser.parse_args()

# Starting
#----------------------------
if prm.prt == True:
    starttime = time.clock()
    print "Starting calculation at",
    print time.strftime("%H:%M:%S on %a %d %b %Y")

# check input
if not os.path.isfile(prm.POSCAR):
    print "\n** ERROR: Initial position file %s was not found." % prm.POSCAR
    sys.exit(0)

if prm.get_super == True:
    # initialize Tmat
    Tmat = np.eye(3)
    # construct Tmat
    if len(map(int, prm.T_mat.split())) == 9:
        Tmat = np.reshape(map(int, prm.T_mat.split()),[3,3])
    elif len(map(int, prm.T_mat.split())) == 3:
        np.fill_diagonal(Tmat,map(int, prm.T_mat.split()))
        # print Tmat
    else:
        print "\n** ERROR: dimenstion of the Transformation mat wrong."
        sys.exit(0)


# if not is_number(prm.symmetry):
#     print "\n** ERROR: symmetry tolerance shoule be a number."
#     sys.exit(0)

# Read information from command line
# First specify location of POSCAR
i_POSCAR=prm.POSCAR.lstrip()

if prm.prt == True:
    print "\nPosition file name: %s " % i_POSCAR
    # print "Symmetry tolerance: %s" % prm.symmetry

# get cell informations
#----------------------------
initial_pos = read(i_POSCAR, format='vasp')
lattice     = initial_pos.get_cell()
positions   = initial_pos.get_scaled_positions()
numbers     = initial_pos.get_atomic_numbers()
# Magnetic moments are not considered in get_spacegroup method
#magmoms = [np.ones(initial_pos.get_number_of_atoms())]
#initial_cell = (lattice, positions, numbers)

if prm.prt == True:
    print '\n==========================================='
    print '\nInitial Structure'
    print "\nLattice Matrix  : (in Angstrom) "
    print lattice
    print "\nAtomic Positions: (in direct coordinate) "
    print positions
    print "\nAtomic numbers  : (for each atom) "
    print numbers
    print '\n==========================================='
# visualize
if prm.visualize == True:
    view(initial_pos)

# get atomic species
fin = [initial_pos.get_chemical_symbols()[0]]
for i in range(len(initial_pos.get_chemical_symbols())):
    if i == len(initial_pos.get_chemical_symbols())-1: break
    if initial_pos.get_chemical_symbols()[i] != initial_pos.get_chemical_symbols()[i+1]:
        fin.append(initial_pos.get_chemical_symbols()[i+1])


# make supercell
#----------------------------
if prm.get_super == True:
    del initial_pos.constraints
    super_cell = make_supercell(initial_pos,Tmat)

    # some post processing...
    #----------------------------
    indices = []
    symb    = []
    for symbol in fin:
        ii = 0
        for atom in super_cell.get_chemical_symbols():
            if atom == symbol:
                indices.append(ii) # get sorted atomic index
                symb.append(symbol) # get sorted atomic symbols
            ii += 1

    pos_new     = [super_cell.get_positions()[i] for i in indices]
    symbls      = symb
    super_cell.set_positions(pos_new)
    super_cell.set_chemical_symbols(symbls)

    lattice     = super_cell.get_cell()
    positions   = super_cell.get_scaled_positions()
    numbers     = super_cell.get_atomic_numbers()

    if prm.prt == True:
        print '\n==========================================='
        print '\nSuper Structure'
        print "\nLattice Matrix  : (in Angstrom) "
        print lattice
        print "\nAtomic Positions: (in direct coordinate) "
        print positions
        print "\nAtomic numbers  : (for each atom) "
        print numbers
        print '\n==========================================='
    # visualize?
    if prm.visualize==True:
        view(super_cell)
    # Write super structure to file
    if prm.over == True:
        write("POSCAR", super_cell, format='vasp')
    else:
        write(i_POSCAR+"_"+"super"+".vasp",super_cell,format='vasp')
    # adding atomic species
    fin = [super_cell.get_chemical_symbols()[0]]
    for i in range(len( super_cell.get_chemical_symbols())):
        if i == len( super_cell.get_chemical_symbols())-1: break
        if  super_cell.get_chemical_symbols()[i] !=  super_cell.get_chemical_symbols()[i+1]:
            fin.append( super_cell.get_chemical_symbols()[i+1])

    if prm.over == True:
        f = open("POSCAR", "r")
    else:
        f = open(i_POSCAR+"_"+"super"+".vasp", "r")
    contents = f.readlines()
    f.close()

    contents.insert(5, ' '.join(fin)+'\n')

    f = open(i_POSCAR+"_"+"super"+".vasp", "w")
    contents = "".join(contents)
    f.write(contents)
    f.close()

    if prm.prt == True:
        if prm.over == False:
            print "\nOutput file name: %s " % str(i_POSCAR+"_"+"super"+".vasp")


# get distance
#----------------------------
if prm.getdist == True:
    atom_list = map(int, prm.atom_list.split())
    # here we use natural: atomic number starting from 1.
    dist = initial_pos.get_distance(atom_list[0]-1,atom_list[1]-1)
    if prm.prt == True:
        sys.stdout.write("\033[1;31m" ) # set color red
        print "\nDistance between atom %d and %d: %f" %(atom_list[0], atom_list[1], dist)
        sys.stdout.write("\033[0;0m") # reset color
    else:
        print dist


if prm.prt == True:
    # Post process
    #----------------------------
    endtime = time.clock()
    runtime = endtime-starttime
    print "\nEnd of calculation."
    print "Program was running for %.2f seconds." % runtime
