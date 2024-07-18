#!/usr/bin/env python
"""
A script which averages a .cube file in three direction to get charge center.
User must specify filename in command line.
eg. python chgcent.py filename.cube "1 2 3" "4 5 6"
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
from ase.io.cube import read_cube_data

# Command line praser
#----------------------------
parser = argparse.ArgumentParser(description='A script to calculate the charge center by integrating along lattice directions.')
parser.add_argument('CHGCAR', nargs='?', default='filename.cube', help="name of the CHGCAR file.")
parser.add_argument('ZVAL', nargs='?', default='1', help="ionic charge of each atom type.")
parser.add_argument('ions_per_type', nargs='?', default='1', help="number of atom of each atom type.")

prm = parser.parse_args()

#-------------------------------------------------------------------------------
np.set_printoptions(formatter={'float': '{: 0.5f}'.format})

starttime = time.time()
print("Starting calculation at",end='')
print(time.strftime("%H:%M:%S on %a %d %b %Y"))


if not os.path.isfile(prm.CHGCAR):
    print("\n** ERROR: Input file %s was not found." % prm.CHGCAR)
    sys.exit(0)

# Read information from command line
# First specify location of CHGCAR
CHGCARfile = prm.CHGCAR.lstrip()


# Open geometry and density class objects
#-----------------------------------------
potl, atoms = read_cube_data(CHGCARfile)
potl = potl*(1.8897259886**3) # convert to e/Angstrom

#sort zval and ions_per_type
#-----------------------------------------
# converting to integer
ZVAL = [int(Z) for Z in ZVAL.split()]
ions_per_type = [int(ions_per_type) for ion in ions_per_type.split()]

assert len(ZVAL) == len(ions_per_type), 'ZVAL and ions_per_type must have the same length.'


print("\nElectronic part")
print("-----------------")
print("Reading file: %s" % CHGCARfile)
#print "Performing average in %s direction" % direction


# Read in lattice parameters and scale factor and atomic positions
#---------------------------------------------
cell = atoms.cell
pos = atoms.get_positions()

# Read Lattice stuff
#---------------------------------------------
dv = np.linalg.det(cell)/(np.array(potl).shape[0]*np.array(potl).shape[1]*np.array(potl).shape[2])
dr = cell/np.array(potl).shape

# Find length of lattice vectors
#--------------------------------
latticelength = np.dot(cell, cell.T).diagonal()
latticelength = latticelength**0.5


# Read in potential data
#------------------------
ngridpts = np.array(potl.shape)
totgridpts = ngridpts.prod()
print("Potential stored on a %dx%dx%d grid" % (ngridpts[0],ngridpts[1],ngridpts[2]))
print("Total number of points is %d" % totgridpts)
print("Reading potential data from file...",end='')
sys.stdout.flush()
print("done.")


# Perform average
#-----------------
average_x = np.zeros(ngridpts[0],float)
average_y = np.zeros(ngridpts[1],float)
average_z = np.zeros(ngridpts[2],float)


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
#  xdiff = latticelength[0]/float(ngridpts[0])
#  ydiff = latticelength[1]/float(ngridpts[1])
#  zdiff = latticelength[2]/float(ngridpts[2])

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

# Total electron numbers
#-------------------
total_elect = np.sum(np.array(potl))*dv
print("\nTotal integrated electron number = %8.5f" % total_elect)
# total_elect = 0
# for i in range(len(ZVAL)):
#     total_elect += int(ions_per_type[i])*int(float(ZVAL[i]))

# Print out average
#-------------------
print("\nCell matrix:")
print(cell)

print("\nCrystal axis length (a,b,c)      = %8.5f , %8.5f , %8.5f" % (latticelength[0],
                                                                      latticelength[1],
                                                                      latticelength[2]))
print("Charge center (crystal axis)     = %8.5f , %8.5f , %8.5f" % (avg_x, avg_y, avg_z))
print("Charge center (Cartisen axis)    = %8.5f , %8.5f , %8.5f" % (avg_x*cell[0,0]+avg_y*cell[1,0]+avg_z*cell[2,0],
                                                                    avg_x*cell[0,1]+avg_y*cell[1,1]+avg_z*cell[2,1],
                                                                    avg_x*cell[0,2]+avg_y*cell[1,2]+avg_z*cell[2,2]))

sys.stdout.write("\033[0;32m") # set color green
print("\nTotal elect dipole moments (e*A) =  %8.5f  %8.5f  %8.5f" % (total_elect*(avg_x*cell[0,0]+avg_y*cell[1,0]+avg_z*cell[2,0]),
                                                                     total_elect*(avg_x*cell[0,1]+avg_y*cell[1,1]+avg_z*cell[2,1]),
                                                                     total_elect*(avg_x*cell[0,2]+avg_y*cell[1,2]+avg_z*cell[2,2])))
sys.stdout.write("\033[0;0m") # reset color

# get ion types and numbers and val
#-------------------
print("\nIonic part")
print("-----------------")
#  print("Reading file: %s" % OUTCARfile)
print("ZVAL          = ", ZVAL)
print("ions per type = ", ions_per_type)


print("\nIonic positions (Cartisen) for each type:")
for i in range(len(ZVAL)):
    print("Type : %i" % int(i))
    print(pos[sum(int(n) for n in (ions_per_type[:i])):sum(int(n) for n in (ions_per_type[:i])) + int(ions_per_type[i])])

# calculate total ionic dipole moments
total_ion_dipole = [0,0,0]
for i in range(len(ZVAL)):
    #print "type : %i" % int(i)
    #print int(float(ZVAL[i]))*pos[sum(int(n) for n in (ions_per_type[:i])):sum(int(n) for n in (ions_per_type[:i])) + int(ions_per_type[i])]
    total_ion_dipole += sum(int(float(ZVAL[i]))*pos[sum(int(n) for n in (ions_per_type[:i])):sum(int(n) for n in (ions_per_type[:i])) + int(ions_per_type[i])])

sys.stdout.write("\033[0;32m") # set color green
print("\nTotal ionic dipole moments (e*A) =  %8.5f  %8.5f  %8.5f" % (total_ion_dipole[0], total_ion_dipole[1], total_ion_dipole[2]))
sys.stdout.write("\033[0;0m") # reset color

print('''
! Total dipole moment comes from:
!     \D_{ion} - \D_{elect}
! So the direction is correct.
! To use C/A unit, simply dot abs(e).''')

sys.stdout.write("\033[1;31m" ) # set color red
print("\nTotal dipole moments (e*A)       =  %8.5f  %8.5f  %8.5f" % (total_ion_dipole[0]-total_elect*(avg_x*cell[0,0]+avg_y*cell[1,0]+avg_z*cell[2,0]),
                                                                     total_ion_dipole[1]-total_elect*(avg_x*cell[0,1]+avg_y*cell[1,1]+avg_z*cell[2,1]),
                                                                     total_ion_dipole[2]-total_elect*(avg_x*cell[0,2]+avg_y*cell[1,2]+avg_z*cell[2,2])))
sys.stdout.write("\033[0;0m") # reset color


# Post process
#-------------------
endtime = time.time()
runtime = endtime-starttime
print("\nEnd of calculation.")
print("Program was running for %.2f seconds." % runtime)
