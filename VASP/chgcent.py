#!/usr/bin/env python
"""
A script which averages a CHGCAR  file in three direction to get charge center.
User must specify filename in command line.
eg. python chgcent.py CHGCAR
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


if not os.path.isfile(sys.argv[1]):
    print "\n** ERROR: Input file %s was not found." % sys.argv[1]
    sys.exit(0)


# Read information from command line
# First specify location of LOCPOT
LOCPOTfile = sys.argv[1].lstrip()


# Open geometry and density class objects
#-----------------------------------------
vasp_charge = VaspChargeDensity(filename = LOCPOTfile)
potl = vasp_charge.chg[-1]
atoms = vasp_charge.atoms[-1]
del vasp_charge


print "\nReading file: %s" % LOCPOTfile
#print "Performing average in %s direction" % direction


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
average_x = np.zeros(ngridpts[0],np.float)
average_y = np.zeros(ngridpts[1],np.float)
average_z = np.zeros(ngridpts[2],np.float)


# At each point, sum over other two indices
for ipt in range(ngridpts[0]):
        average_x[ipt] = potl[ipt,:,:].sum()


for ipt in range(ngridpts[1]):
        average_y[ipt] = potl[:,ipt,:].sum()


for ipt in range(ngridpts[2]):
        average_z[ipt] = potl[:,:,ipt].sum()


# Scale by size of area element in the plane,
# gives unit e/Ang. I.e. integrating the resulting
# CHG_dir file should give the total charge.
area_x = np.linalg.det([ (cell[1,1], cell[1,2] ),(cell[2,1], cell[2,2])])
area_y = np.linalg.det([ (cell[0,0], cell[0,2] ),(cell[2,0], cell[2,2])])
area_z = np.linalg.det([ (cell[0,0], cell[0,1] ),(cell[1,0], cell[1,1])])


dA_x = area_x/(ngridpts[1]*ngridpts[2])
dA_y = area_y/(ngridpts[0]*ngridpts[2])
dA_z = area_z/(ngridpts[0]*ngridpts[1])


average_x *= dA_x
average_y *= dA_y
average_z *= dA_z


# get unit length using NGXF NGYF and NGZF
# position start at 0
xdiff = latticelength[0]/float(ngridpts[0])
ydiff = latticelength[1]/float(ngridpts[1])
zdiff = latticelength[2]/float(ngridpts[2])

avg_x=avg_y=avg_z=0.0
# Calculate averaged position
for i in range(ngridpts[0]):
    x = float(i)/ngridpts[0]
    avg_x += x*average_x[i]/average_x.sum()


for i in range(ngridpts[1]):
    y = float(i)/ngridpts[1]
    avg_y += y*average_y[i]/average_y.sum()


for i in range(ngridpts[2]):
    z = float(i)/ngridpts[2]
    avg_z += z*average_z[i]/average_z.sum()


# Print out average
#-------------------
print "\nCell matrix:"
print cell

print "\nCrystal axis length (a,b,c) = %15.8g , %15.8g , %15.8g" % (latticelength[0],latticelength[1], latticelength[2])

print "\nCharge center (crystal axis) = %15.8g , %15.8g , %15.8g" % (avg_x, avg_y, avg_z)

print "\nCharge center (Cartisen axis) = %15.8g , %15.8g , %15.8g" % (avg_x*cell[0,0] + avg_y*cell[1,0] + avg_z*cell[2,0], avg_x*cell[0,1] + avg_y*cell[1,1] + avg_z*cell[2,1], avg_x*cell[0,2] + avg_y*cell[1,2] + avg_z*cell[2,2])

# Post process
#-------------------
endtime = time.clock()
runtime = endtime-starttime
print "\nEnd of calculation."
print "Program was running for %.2f seconds." % runtime

