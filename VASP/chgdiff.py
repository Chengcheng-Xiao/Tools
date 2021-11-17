#!/usr/bin/env python
"""
A script which reads CHGCAR files and takes the difference
between the first density and the following ones.
Output is to a file called CHGDIFF.
Depends on ase
"""
from __future__ import print_function
import argparse
import os
import sys
import numpy
import time
from ase.calculators.vasp import VaspChargeDensity

# Command line praser
#----------------------------
parser = argparse.ArgumentParser(description='A script to calculate the differential charge density.')
parser.add_argument('CHGCAR', nargs='*', help="name of the CHGCAR files.")
# parser.add_argument('OUTCAR', nargs='?', default='OUTCAR', help="name of the OUTCAR file.")

prm = parser.parse_args()


# Find out how many arguments were on the command line,
# all but the first two should contain files with densities
# to subtract
# nsubtract = len(sys.argv)-2
nsubtract = len(prm.CHGCAR)
if not nsubtract >= 1:
    print("\n** ERROR: Must specify name of at least two files on command line.")
    print("eg. chgdiff.py CHGCAR1 CHGCAR2  [CHGCAR3 ...]")
    print("The reference density is taken from the first filename.")
    print("The densities in the files after this will be subtracted from the reference.")
    sys.exit(0)

# Check that files exist
for name in prm.CHGCAR:
    if not os.path.isfile(name):
        print("\n** ERROR: Input file %s was not found." % name)
        sys.exit(0)

starttime = time.time()
print("Starting calculation at",end='')
print(time.strftime("%H:%M:%S on %a %d %b %Y \n"))



# Read information from command line
# First specify location of CHGCAR file with reference density
CHGCARfile1 = sys.argv[1].lstrip()

# Open geometry and density class objects
#-----------------------------------------
print("Reading data from file %s ..." % CHGCARfile1, end='')
sys.stdout.flush()
vasp_charge1 = VaspChargeDensity(filename = CHGCARfile1)
chg1 = vasp_charge1.chg[-1]
atoms1 = vasp_charge1.atoms[-1]
del vasp_charge1
print("done.")

chgdiff=chg1
for CHGCARfile2 in sys.argv[2:]:
    CHGCARfile2 = CHGCARfile2.strip()
    print("Reading data from file %s ..." % CHGCARfile2, end='')
    sys.stdout.flush()
    vasp_charge2 = VaspChargeDensity(filename = CHGCARfile2)
    chg2 = vasp_charge2.chg[-1]
    del vasp_charge2
    print("done.")

    # Make sure that the second data set is on the same grid
    #--------------------------------------------------------
    if chg2.shape != chg1.shape:
       print("\n**ERROR: Two sets of data are not on the same grid.")
       print("Data from file %s on %dx%dx%d grid." % (CHGCARfile1,chg1.shape[0],chg1.shape[1],chg1.shape[2]))
       print("Data from file %s on %dx%dx%d grid.\n" % (CHGCARfile2,chg2.shape[0],chg2.shape[1],chg2.shape[2]))
       sys.exit(0)
    else:
       print("Subtracting %s from file %s ..." % (CHGCARfile2, CHGCARfile1), end='')
       sys.stdout.flush()

    # Take difference
    #-----------------
    chgdiff=chgdiff-chg2
    print("done.")


vasp_charge_diff = VaspChargeDensity(filename=None)
vasp_charge_diff.atoms=[atoms1,]
vasp_charge_diff.chg=[chgdiff,]

# Print out charge density
#--------------------------
# Check whether CHGDIFF exists
if os.path.isfile("./CHGDIFF.vasp"):
    print("\n**WARNING: A file called CHGDIFF already exists in this directory.")
    if float(sys.version.split()[0][:3]) < 3.0:
        yesno=raw_input("Type y to continue and overwrite it, any other key to stop\n")
    else:
        yesno=input("Type y to continue and overwrite it, any other key to stop\n")
    if yesno!="y":
      sys.exit(0)


print("Writing density difference data to file CHGDIFF ...", end='')
sys.stdout.flush()
vasp_charge_diff.write(filename="CHGDIFF.vasp",format="chgcar")

if float(sys.version.split()[0][:3]) < 3.0:
    # adding atomic species
    fin = [atoms1.get_chemical_symbols()[0]]
    for i in range(len(atoms1.get_chemical_symbols())):
        if i == len(atoms1.get_chemical_symbols())-1: break
        if atoms1.get_chemical_symbols()[i] != atoms1.get_chemical_symbols()[i+1]:
            fin.append(atoms1.get_chemical_symbols()[i+1])

    for chgfiles in ['CHGDIFF.vasp']:
        f = open(chgfiles, "r")
        contents = f.readlines()
        f.close()

        contents.insert(5, ' '.join(fin)+'\n')

        f = open(chgfiles, "w")
        contents = "".join(contents)
        f.write(contents)
        f.close()

print("done.")

endtime = time.time()
runtime = endtime-starttime
print("\nEnd of calculation.")
print("Program was running for %.2f seconds." % runtime)
