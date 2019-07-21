#!/usr/bin/env python
"""
A script to find the the symmetry of unitcell and return refined structures.
Usage:
./impose_sym.py POSCAR 1E-5

Depends on ase and spglib
"""

import spglib
from ase.io import read, write
import numpy as np
import os
import sys
import datetime
import time

# is digit?
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


# Starting
#----------------------------
starttime = time.clock()
print "Starting calculation at",
print time.strftime("%H:%M:%S on %a %d %b %Y")

if not os.path.isfile(sys.argv[1]):
    print "\n** ERROR: Initial position file %s was not found." % sys.argv[1]
    sys.exit(0)

if not is_number(sys.argv[2]):
    print "\n** ERROR: symmetry tolerance shoule be a number."
    sys.exit(0)

# Read information from command line
# First specify location of POSCAR
i_POSCAR=sys.argv[1].lstrip()

print "\nPosition file name: %s " % i_POSCAR
print "Symmetry tolerance: %s" % sys.argv[2]


# get cell informations
#----------------------------
initial_pos = read(i_POSCAR)
lattice     = initial_pos.get_cell()
positions   = initial_pos.get_scaled_positions()
numbers     = initial_pos.get_atomic_numbers()
# Magnetic moments are not considered in get_spacegroup method
#magmoms = [np.ones(initial_pos.get_number_of_atoms())]
initial_cell = (lattice, positions, numbers)

# find symmetry
#----------------------------
spacegroup = spglib.get_spacegroup(initial_cell, symprec=float(sys.argv[2]))
print "Spacegroup: %s" % spacegroup

f_lattice, f_positions, f_numbers= spglib.refine_cell(initial_cell, symprec=float(sys.argv[2]))
print "\nLattice Matrix  : (in Angstrom) "
print f_lattice
print "\nAtomic_positions: (in direct coordinate) "
print f_positions
print "\nAtomic numbers  : (for each atom) "
print f_numbers


# set cell informations
#----------------------------
final_pos = read(i_POSCAR)
final_pos.set_cell(f_lattice)
final_pos.set_scaled_positions(f_positions)
final_pos.set_atomic_numbers(f_numbers)

write(i_POSCAR+"_"+sys.argv[2]+".vasp",final_pos,format='vasp')

print "\nOutput file name: %s " % str(i_POSCAR+"_"+sys.argv[2]+".vasp")
# Post process
#----------------------------
endtime = time.clock()
runtime = endtime-starttime
print "\nEnd of calculation."
print "Program was running for %.2f seconds." % runtime
