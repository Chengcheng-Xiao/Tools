#!/usr/bin/env python
"""
A script to read POSCAR/CONTCAR files and visualize them using ase.gui
Depends on ase
"""

from __future__ import print_function
import sys, os
import argparse
from ase.io import read
from ase.visualize import view

# Command line praser
#----------------------------
parser = argparse.ArgumentParser(description='A script to view atomic structure.')
parser.add_argument('POSCAR', nargs='*', help="name of the CHGCAR files.")

prm = parser.parse_args()

# Find out how many arguments were on the command line,
# nsubtract = len(sys.argv)-2
nsubtract = len(prm.POSCAR)
if nsubtract > 1:
    print("\n** ERROR: Can only view one file at a time.")
    print("eg. view_atoms.py POSCAR")
    sys.exit(0)

# Check that files exist
for name in prm.POSCAR:
    if not os.path.isfile(name):
        print("\n** ERROR: Input file %s was not found." % name)
        sys.exit(0)

# Read information from command line
File = prm.POSCAR[0].lstrip()

atoms = read(File, format='vasp')
view(atoms)
