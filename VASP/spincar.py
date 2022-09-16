#!/usr/bin/env python
"""
A script to extract spindensity form CHGCAR file.
User must specify filename in command line.
eg. python spincar.py CHGCAR
** Does NOT work with Molecular dynamics output.
Depends on ase
"""

from __future__ import print_function
import argparse
import os
import sys
import numpy as np
import math
import string
import datetime
import time
from ase.calculators.vasp import VaspChargeDensity


# Command line praser
#----------------------------
parser = argparse.ArgumentParser(description='A script to calculate the spin density.')
parser.add_argument('CHGCAR', nargs='*', help="name of the CHGCAR files.")

prm = parser.parse_args()

starttime = time.time()
print("Starting calculation at",end='')
print(time.strftime("%H:%M:%S on %a %d %b %Y"))


# Find out how many arguments were on the command line,
nsubtract = len(prm.CHGCAR)
if not nsubtract == 1:
    print("\n** ERROR: Must specify the name of file on command line.")
    print("eg. chgdiff.py CHGCAR")
    print("Only one file name are allowed")
    sys.exit(0)

if not os.path.isfile(prm.CHGCAR[0]):
    print("\n** ERROR: Input file %s was not found." % prm.CHGCAR[0])
    sys.exit(0)


# Read information from command line
# First specify location of CHGCAR
CHGCARfile = prm.CHGCAR[0].lstrip()


# Open geometry and density class objects
#-----------------------------------------
vasp_charge = VaspChargeDensity(filename = CHGCARfile)
if len(vasp_charge.chgdiff) == 3:
    spin=3
    print("\nReading spin orbital potential data from file %s ... " % CHGCARfile,end='')
    sys.stdout.flush()
    atoms = vasp_charge.atoms[-1]
    potl_total = vasp_charge.chg[-1]
    potl_x = vasp_charge.chgdiff[0]
    potl_y = vasp_charge.chgdiff[1]
    potl_z = vasp_charge.chgdiff[2]
    if not potl_total.shape == potl_x.shape == potl_y.shape == potl_z.shape:
        print("\n**ERROR: Two sets of data are not on the same grid.")
        print("Data from data block 1 on %dx%dx%d grid." % (potl_total.shape[0],potl_total.shape[1],potl_total.shape[2]))
        print("Data from data block 2 on %dx%dx%d grid." % (potl_x.shape[0],potl_x.shape[1],potl_x.shape[2]))
        print("Data from data block 3 on %dx%dx%d grid." % (potl_y.shape[0],potl_y.shape[1],potl_y.shape[2]))
        print("Data from data block 4 on %dx%dx%d grid." % (potl_z.shape[0],potl_z.shape[1],potl_z.shape[2]))
        sys.exit(0)
    else:
        print('done.')

elif len(vasp_charge.chgdiff) == 1:
    spin=2
    print("\nReading spin polarized potential data from file %s ... " % CHGCARfile,end='')
    sys.stdout.flush()
    atoms = vasp_charge.atoms[-1]
    potl_updown = vasp_charge.chg[-1]
    potl_upmdown = vasp_charge.chgdiff[0]
    if not potl_updown.shape == potl_upmdown.shape:
        print("\n**ERROR: Two sets of data are not on the same grid.")
        print("Data from data block 1 on %dx%dx%d grid." % (potl_updown.shape[0],potl_updown.shape[1],potl_updown.shape[2]))
        print("Data from data block 2 on %dx%dx%d grid." % (potl_updown.shape[0],potl_upmdown.shape[1],potl_upmdown.shape[2]))
        sys.exit(0)
    else:
        # get density for each spin channels
        potl_up = (potl_updown+potl_upmdown)/2
        potl_down = (potl_updown-potl_upmdown)/2
        print('done.')

elif len(vasp_charge.chgdiff) == 0:
    spin=1
    print("\nReading spin paired potential data from file %s ... " % CHGCARfile, end='')
    sys.stdout.flush()
    print("\nFile only contains one block of data, check if is polarized calculation.")
    sys.exit(0)

else:
    print("\n** ERROR: CHGCAR data block error, total of %s data blocks was found." % len(vasp_charge.chg))
    sys.exit(0)

del vasp_charge


# Exctract data
#------------------------
if spin == 3: # SOC
    chg_total = VaspChargeDensity(filename=None)
    chg_total.atoms=[atoms,]
    chg_total.chg=[potl_total,]

    chg_x = VaspChargeDensity(filename=None)
    chg_x.atoms=[atoms,]
    chg_x.chg=[potl_x,]

    chg_y = VaspChargeDensity(filename=None)
    chg_y.atoms=[atoms,]
    chg_y.chg=[potl_y,]

    chg_z = VaspChargeDensity(filename=None)
    chg_z.atoms=[atoms,]
    chg_z.chg=[potl_z,]

    # Print out charge density
    #--------------------------
    # Check whether CHGDIFF exists
    for chgfiles in ['SPINCAR_total.vasp','SPINCAR_x.vasp','SPINCAR_y.vasp','SPINCAR_z.vasp']:
        if os.path.isfile("./"+i):
            print("\n**WARNING: A file called CHGDIFF already exists in this directory.")
            if float(sys.version.split()[0][:3]) < 3.0:
                yesno=raw_input("Type y to continue and overwrite it, any other key to stop\n")
            else:
                yesno=input("Type y to continue and overwrite it, any other key to stop\n")
            if yesno!="y":
                sys.exit(0)
    # Writing data ...
    print("Writing density difference data to file ...", end='')
    sys.stdout.flush()
    chg_total.write(filename="SPINCAR_total.vasp",format="chgcar")
    chg_x.write(filename="SPINCAR_x.vasp",format="chgcar")
    chg_y.write(filename="SPINCAR_y.vasp",format="chgcar")
    chg_z.write(filename="SPINCAR_z.vasp",format="chgcar")
    # adding atomic species
    fin = [atoms.get_chemical_symbols()[0]]
    for i in range(len(atoms.get_chemical_symbols())):
        if i == len(atoms.get_chemical_symbols())-1: break
        if atoms.get_chemical_symbols()[i] != atoms.get_chemical_symbols()[i+1]:
            fin.append(atoms.get_chemical_symbols()[i+1])

    for chgfiles in ['SPINCAR_total.vasp','SPINCAR_x.vasp','SPINCAR_y.vasp','SPINCAR_z.vasp']:
        f = open(chgfiles, "r")
        contents = f.readlines()
        f.close()

        contents.insert(5, ' '.join(fin)+'\n')

        f = open(chgfiles, "w")
        contents = "".join(contents)
        f.write(contents)
        f.close()
    print("done.")

if spin == 2: # Spin polarized
    chg_updown = VaspChargeDensity(filename=None)
    chg_updown.atoms=[atoms,]
    chg_updown.chg=[potl_updown,]

    chg_upmdown = VaspChargeDensity(filename=None)
    chg_upmdown.atoms=[atoms,]
    chg_upmdown.chg=[potl_upmdown,]

    chg_up = VaspChargeDensity(filename=None)
    chg_up.atoms=[atoms,]
    chg_up.chg=[potl_up,]

    chg_down = VaspChargeDensity(filename=None)
    chg_down.atoms=[atoms,]
    chg_down.chg=[potl_down,]

    # Print out charge density
    #--------------------------
    # Check whether CHGDIFF exists
    for chgfiles in ['SPINCAR_up+down.vasp','SPINCAR_up-down.vasp','SPINCAR_up.vasp','SPINCAR_down.vasp']:
        if os.path.isfile("./"+chgfiles):
            print("\n**WARNING: A file called ' " + chgfiles + " ' already exists in this directory.")
            if float(sys.version.split()[0][:3]) < 3.0:
                yesno=raw_input("Type y to continue and overwrite it, any other key to stop\n")
            else:
                yesno=input("Type y to continue and overwrite it, any other key to stop\n")
            if yesno!="y":
                sys.exit(0)
    # Writing data ...
    print("\nWriting density difference data to file ...", end='')
    sys.stdout.flush()
    chg_updown.write(filename="SPINCAR_up+down.vasp",format="chgcar")
    chg_upmdown.write(filename="SPINCAR_up-down.vasp",format="chgcar")
    chg_up.write(filename="SPINCAR_up.vasp",format="chgcar")
    chg_down.write(filename="SPINCAR_down.vasp",format="chgcar")

    if float(sys.version.split()[0][:3]) < 3.0:
        # adding atomic species
        fin = [atoms.get_chemical_symbols()[0]]
        for i in range(len(atoms.get_chemical_symbols())):
            if i == len(atoms.get_chemical_symbols())-1: break
            if atoms.get_chemical_symbols()[i] != atoms.get_chemical_symbols()[i+1]:
                fin.append(atoms.get_chemical_symbols()[i+1])

        for chgfiles in ['SPINCAR_up+down.vasp','SPINCAR_up-down.vasp','SPINCAR_up.vasp','SPINCAR_down.vasp']:
            f = open(chgfiles, "r")
            contents = f.readlines()
            f.close()

            contents.insert(5, ' '.join(fin)+'\n')

            f = open(chgfiles, "w")
            contents = "".join(contents)
            f.write(contents)
            f.close()
    print("done.")

# Post process
#-------------------
endtime = time.time()
runtime = endtime-starttime
print("\nEnd of calculation.")
print("Program was running for %.2f seconds." % runtime)
