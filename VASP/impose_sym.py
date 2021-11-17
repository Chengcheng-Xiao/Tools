#!/usr/bin/env python
"""
A script to find the the symmetry of unitcell and return refined structures.

Depends on ase and spglib
"""

from __future__ import print_function
import argparse
import spglib
from ase.io import read, write
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
# Do we want to convert the final structure into use primitive cell?
parser = argparse.ArgumentParser(description='A script to find the the symmetry of unitcell and return refined structures.')
parser.add_argument('-noi','--no_ideal', action="store_true", default=False, dest='no_ideal',
                    help='Do not impose symmetry. If set, only find primitive cell without change atomic positions. (Default=False)')
parser.add_argument('-v','--visualize', action="store_true", default=False, dest="visualize",
                    help='Use ASE-GUI to visualize structure (Default=False)')
parser.add_argument('-s','--symmetry', action="store", default=float(1E-5), dest="symmetry",
                    help='symmetry tolerance for finding symmetry (Default=1E-5)')
parser.add_argument('-c','--POSCAR', action="store", default=str('POSCAR'), dest="POSCAR",
                    help='input file name (Default=POSCAR)')
prm = parser.parse_args()

# Starting
#----------------------------
starttime = time.time()
print("Starting calculation at", end='')
print(time.strftime("%H:%M:%S on %a %d %b %Y"))

# check input
if not os.path.isfile(prm.POSCAR):
    print("\n** ERROR: Initial position file %s was not found." % prm.POSCAR)
    sys.exit(0)

if not is_number(prm.symmetry):
    print("\n** ERROR: symmetry tolerance shoule be a number.")
    sys.exit(0)

# Read information from command line
# First specify location of POSCAR
i_POSCAR=prm.POSCAR.lstrip()

print("\nPosition file name: %s " % i_POSCAR)
print("Symmetry tolerance: %s" % prm.symmetry)

# get cell informations
#----------------------------
initial_pos = read(i_POSCAR)
lattice     = initial_pos.get_cell()
positions   = initial_pos.get_scaled_positions()
numbers     = initial_pos.get_atomic_numbers()
# Magnetic moments are not considered in get_spacegroup method
#magmoms = [np.ones(initial_pos.get_number_of_atoms())]
initial_cell = (lattice, positions, numbers)

print('\n===========================================')
print('\nInitial Structure')
print("\nLattice Matrix  : (in Angstrom) ")
print(lattice)
print("\nAtomic Positions: (in direct coordinate) ")
print(positions)
print("\nAtomic numbers  : (for each atom) ")
print(numbers)
print('\n===========================================')

# visualize
if prm.visualize==True:
    view(initial_pos)

# find symmetry
#----------------------------
spacegroup = spglib.get_spacegroup(initial_cell, symprec=float(prm.symmetry))
print("Spacegroup: %s" % spacegroup)

print('\n===========================================')
print("\nFinding primitive cell...")

# impoer or not?
if prm.no_ideal==False:
    print('\nImposing symmetry...')
    f_lattice, f_positions, f_numbers = spglib.standardize_cell(initial_cell, to_primitive=True, no_idealize=False, symprec=float(prm.symmetry))
else:
    print('\nsymmetry not imposed')
    f_lattice, f_positions, f_numbers = spglib.standardize_cell(initial_cell, to_primitive=True, no_idealize=True, symprec=float(prm.symmetry))

# Final structure
print("\nLattice Matrix  : (in Angstrom) ")
print(f_lattice)
print("\nAtomic Positions: (in direct coordinate) ")
print(f_positions)
print("\nAtomic numbers  : (for each atom) ")
print(f_numbers)
print('\n===========================================')


# set cell informations
#----------------------------
final_pos = Atoms(numbers=f_numbers,
                  pbc=True,
                  cell=f_lattice,
                  scaled_positions=f_positions)
# visualize
if prm.visualize==True:
    view(final_pos)

# Write final structure to file
write(i_POSCAR+"_"+str(prm.symmetry)+".vasp",final_pos,format='vasp')

if float(sys.version.split()[0][:3]) < 3.0:
    # adding atomic species
    fin = [final_pos.get_chemical_symbols()[0]]
    for i in range(len(final_pos.get_chemical_symbols())):
        if i == len(final_pos.get_chemical_symbols())-1: break
        if final_pos.get_chemical_symbols()[i] != final_pos.get_chemical_symbols()[i+1]:
            fin.append(final_pos.get_chemical_symbols()[i+1])

    f = open(i_POSCAR+"_"+str(prm.symmetry)+".vasp", "r")
    contents = f.readlines()
    f.close()

    contents.insert(5, ' '.join(fin)+'\n')

    f = open(i_POSCAR+"_"+str(prm.symmetry)+".vasp", "w")
    contents = "".join(contents)
    f.write(contents)
    f.close()

print("\nOutput file name: %s " % str(i_POSCAR+"_"+str(prm.symmetry)+".vasp"))
# Post process
#----------------------------
endtime = time.time()
runtime = endtime-starttime
print("\nEnd of calculation.")
print("Program was running for %.2f seconds." % runtime)
