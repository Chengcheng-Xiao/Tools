#!/usr/bin/env python
"""
A script to read POSCAR/CONTCAR files and visualize them using ase.gui
Depends on ase
"""
import sys, os
from ase.io import read
from ase.visualize import view

# Find out how many arguments were on the command line,
nsubtract = len(sys.argv)-2
if nsubtract >= 1:
    print "\n** ERROR: Only support one file."
    print "eg. view_atoms.py POSCAR"
    sys.exit(0)

# Check that files exist
for name in sys.argv[1:]:
    if not os.path.isfile(name):
        print "\n** ERROR: Input file %s was not found." % name
        sys.exit(0)

# Read information from command line
File = sys.argv[1].lstrip()

atoms = read(File, format='vasp')
view(atoms)
