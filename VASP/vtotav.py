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
import argparse
from ase.calculators.vasp import VaspChargeDensity
from scipy import interpolate

# Command line praser
#----------------------------
parser = argparse.ArgumentParser(description='A script to calculate planar and macroscopic average along certain axis.')
parser.add_argument('-c', action="store", default="LOCPOT", dest="LOCPOT",
                    help='Input file name.')
parser.add_argument('-chg','--chg', action="store_true", default=False, dest="chg",
                    help='Is density quantity? (Default=False)')
parser.add_argument('-d', action="store", default="z", dest="dir",
                    help='direction(Default=Z)')
parser.add_argument('-macro', action="store_true", default=False, dest="macro",
                    help='macroscopic average? (Default=False). \
                    Will print out intrpolated planar and macro average. ')
parser.add_argument('--len', action="store", default="0.1", type=np.float, dest="macro_len",
                    help='atomic plane length. (in $\AA$)')

prm = parser.parse_args()


starttime = time.clock()
print "Starting calculation at",
print time.strftime("%H:%M:%S on %a %d %b %Y")

# if len(sys.argv) != 3:
#     print "\n** ERROR: Must specify name of file and direction on command line."
#     print "eg. vtotav.py LOCPOT z."
#     sys.exit(0)

if not os.path.isfile(prm.LOCPOT):
    print "\n** ERROR: Input file %s was not found." % prm.LOCPOT
    sys.exit(0)

# Read information from command line
# First specify location of LOCPOT
LOCPOTfile = prm.LOCPOT.lstrip()

# Next the direction to make average in
# input should be x y z, or X Y Z. Default is Z.
allowed = "xyzXYZ"
direction = prm.dir.lstrip()
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

# if 'LOCPOT' in LOCPOTfile:
if not prm.chg:
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

if prm.macro == True:
    average_m = np.zeros(ngridpts[idir]*50,np.float)
    interval = prm.macro_len
    dr = np.sqrt((cell[idir]**2).sum())/ngridpts[idir]
    dr_interp = dr/50.
    number = np.floor(np.round(interval/dr)/2)
    number_interp = np.floor(np.round(interval/dr_interp)/2)

    # linear interpolate (50X dense)
    #-------------------
    # x = np.linspace(0, dr*(ngridpts[idir]+1), ngridpts[idir], endpoint=False)
    # y = average
    # xvals = np.linspace(0, dr*(ngridpts[idir]+1), ngridpts[idir]*50, endpoint=False)
    # average_interp = np.interp(xvals, x, y)

    # cubic spline (50X dense)
    #-------------------
    x = np.linspace(0, dr*(ngridpts[idir]+1), ngridpts[idir], endpoint=False)
    y = average
    tck = interpolate.splrep(x, y, s=0)
    xvals = np.linspace(0, dr*(ngridpts[idir]+1), ngridpts[idir]*50, endpoint=False)
    average_interp = interpolate.splev(xvals, tck, der=0)

    # another cubic spline (50X dense)
    #-------------------
    # x = np.linspace(0, dr*(ngridpts[idir]+1), ngridpts[idir], endpoint=False)
    # y = average
    # tck = interpolate.interp1d(x, y, kind='cubic')
    # xvals = np.linspace(0, dr*(ngridpts[idir]), ngridpts[idir]*50, endpoint=False)
    # average_interp = tck(xvals)


    for ipt in range(ngridpts[idir]*50):
        average_m[ipt] = np.array([average_interp[i%(ngridpts[idir]*50)] \
        for i in np.array(np.arange(ipt-number_interp,ipt+number_interp+1,1),dtype=np.int)]).sum()

    # original no interpolation method
    #-------------------
    # for ipt in range(ngridpts[idir]):
    #     average_m[ipt] = np.array([average[i%ngridpts[idir]] \
    #     for i in np.array(np.arange(ipt-number,ipt+number+1,1),dtype=np.int)]).sum()
    #     print [i%ngridpts[idir] \
    #     for i in np.array(np.arange(ipt-number,ipt+number+1,1),dtype=np.int)]

    # # interpolate again?
    #-------------------
    # for ii in range(1):
    #     average_mm = np.zeros(ngridpts[idir]*50,np.float)
    #     for ipt in range(ngridpts[idir]*50):
    #         average_mm[ipt] = np.array([average_m[i%(ngridpts[idir]*50)] \
    #         for i in np.array(np.arange(ipt-number_interp,ipt+number_interp,1),dtype=np.int)]).sum()
    #     average_m = average_mm/(2*number_interp+1)

    average = average_m /(2*number_interp+1)

# Print out average macro
#-------------------
if prm.macro == True:
    averagefile = LOCPOTfile + filesuffix + '_macro'
    print "Writing macroscopic averaged data to file %s..." % averagefile,
    sys.stdout.flush()
    outputfile = open(averagefile,"w")
    if 'LOCPOT' in LOCPOTfile:
        outputfile.write("#  Distance(Ang)     Potential(eV)\n")
    else:
        outputfile.write("#  Distance(Ang)     Chg. density (e/Ang)\n")
    # (50X dense)
    xdiff = latticelength[idir]/float(ngridpts[idir]*50)
    for i in range(ngridpts[idir]*50):
        x = i*xdiff
        outputfile.write("%15.8g %15.8g\n" % (x,average[i]))
    outputfile.close()
    print "done."

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
if prm.macro == True:
    # (50X dense)
    xdiff = latticelength[idir]/float(ngridpts[idir]*50)
    for i in range(ngridpts[idir]*50):
        x = i*xdiff
        outputfile.write("%15.8g %15.8g\n" % (x,average_interp[i]))
else:
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
