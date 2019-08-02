#!/usr/bin/env python
"""
A script which averages a CHGCAR or LOCPOT file in one direction to make a 1D curve.
User must specify filename and direction on command line.
Depends on ase
"""

import os
import sys
import numpy as np
import math
import string
import datetime
import time
from ase.calculators.vasp import VaspChargeDensity

starttime = time.clock()
print "Starting calculation at",
print time.strftime("%H:%M:%S on %a %d %b %Y")

if len(sys.argv) != 3:
    print "\n** ERROR: Must specify name of file and direction on command line."
    print "eg. vtotav.py LOCPOT z."
    sys.exit(0)

if not os.path.isfile(sys.argv[1]):
    print "\n** ERROR: Input file %s was not found." % sys.argv[1]
    sys.exit(0)

# Read information from command line
# First specify location of LOCPOT
LOCPOTfile = sys.argv[1].lstrip()

# Next the direction to make average in
# input should be x y z, or X Y Z. Default is Z.
allowed = "xyzXYZ"
direction = sys.argv[2].lstrip()
if allowed.find(direction) == -1 or len(direction)!=1 :
    print "** WARNING: The direction was input incorrectly."
    print "** Setting to z-direction by default."
if direction.islower():
    direction = direction.upper()
filesuffix = "_%s" % direction

# Open geometry and density class objects
#-----------------------------------------
vasp_charge = VaspChargeDensity(filename = LOCPOTfile)
potl = vasp_charge.chg[-1]
atoms = vasp_charge.atoms[-1]
del vasp_charge

# For LOCPOT files we multiply by the volume to get back to eV
if 'LOCPOT' in LOCPOTfile:
    potl=potl*atoms.get_volume()

print "\nReading file: %s" % LOCPOTfile
print "Performing average in %s direction" % direction

# Read in lattice parameters and scale factor
#---------------------------------------------
cell = atoms.cell

# Find length of lattice vectors
#--------------------------------
latticelength = np.dot(cell, cell.T).diagonal()
latticelength = latticelength**0.5

# Read in potential data
#------------------------
ngridpts = np.array(potl.shape)
totgridpts = ngridpts.prod()
print "Potential stored on a %dx%dx%d grid" % (ngridpts[0],ngridpts[1],ngridpts[2])
print "Total number of points is %d" % totgridpts
print "Reading potential data from file...",
sys.stdout.flush()
print "done."

# Perform average
#-----------------
if direction=="X":
    idir = 0
    a = 1
    b = 2
elif direction=="Y":
    a = 0
    idir = 1
    b = 2
else:
    a = 0
    b = 1
    idir = 2
a = (idir+1)%3
b = (idir+2)%3
# At each point, sum over other two indices
average = np.zeros(ngridpts[idir],np.float)
for ipt in range(ngridpts[idir]):
    if direction=="X":
        average[ipt] = potl[ipt,:,:].sum()
    elif direction=="Y":
        average[ipt] = potl[:,ipt,:].sum()
    else:
        average[ipt] = potl[:,:,ipt].sum()

if 'LOCPOT' in LOCPOTfile:
    # Scale by number of grid points in the plane.
    # The resulting unit will be eV.
    average /= ngridpts[a]*ngridpts[b]
else:
    # Scale by size of area element in the plane,
    # gives unit e/Ang. I.e. integrating the resulting
    # CHG_dir file should give the total charge.
    area = np.linalg.det([ (cell[a,a], cell[a,b] ),
                           (cell[b,a], cell[b,b])])
    dA = area/(ngridpts[a]*ngridpts[b])
    average *= dA

# Print out average
#-------------------
averagefile = LOCPOTfile + filesuffix
print "Writing averaged data to file %s..." % averagefile,
sys.stdout.flush()
outputfile = open(averagefile,"w")
if 'LOCPOT' in LOCPOTfile:
    outputfile.write("#  Distance(Ang)     Potential(eV)\n")
else:
    outputfile.write("#  Distance(Ang)     Chg. density (e/Ang)\n")
xdiff = latticelength[idir]/float(ngridpts[idir])
for i in range(ngridpts[idir]):
    x = i*xdiff
    outputfile.write("%15.8g %15.8g\n" % (x,average[i]))
outputfile.close()
print "done."

endtime = time.clock()
runtime = endtime-starttime
print "\nEnd of calculation."
print "Program was running for %.2f seconds." % runtime
