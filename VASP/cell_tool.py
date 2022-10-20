#!/usr/bin/env python
"""
A script to make supercell and find atomic distance
"""

from __future__ import print_function
import argparse
from ase.io import read, write
from ase.build.supercells import make_supercell
from ase import Atoms
from ase.visualize import view
import numpy as np
import os
import sys
import datetime
import time

# is digit? function
def check_int(s):
    if s[0] in ('-', '+'):
        return s[1:].isdigit()
    return s.isdigit()

# Command line praser
#----------------------------
parser = argparse.ArgumentParser(description='A script to make supercell(transformed supercell) and measure atomic distances.')
parser.add_argument('-o','--overwrite', action="store_true", default=False, dest='over',
                    help='Overwrite POSCAR?. (Default=False)')
parser.add_argument('-vi','--visualize', action="store_true", default=False, dest="visualize",
                    help='Use ASE-GUI to visualize structure (Default=False)')
parser.add_argument('-v', '--verbose', action="store_true", default=False, dest="prt",
                    help='print out info? (Default=False)')
parser.add_argument('-c','--POSCAR', action="store", default=str('POSCAR'), dest="POSCAR",
                    help='Input file name (Default=POSCAR)')
parser.add_argument('-d', action="store_true", default=False, dest="get_super",
                    help='Make supercell? (Default=False)')
parser.add_argument('--dim', action="store", default="1 0 0 0 1 0 0 0 1", dest="T_mat",
                    help='Transformation matrix (full matrix"1,0,0,0,1,0,0,0,1" or diagonal terms "1 1 1")')
parser.add_argument('-di', action="store_true", default=False, dest="getdist",
                    help='Get distance? (Default=False)')
parser.add_argument('--dis', action="store", default="1 2", dest="atom_list",
                    help='Distance between atoms a and b, starting from 1 (Default="1 2")')

prm = parser.parse_args()

# Starting
#----------------------------
if prm.prt == True:
    starttime = time.time()
    print("Starting calculation at", end='')
    print(time.strftime("%H:%M:%S on %a %d %b %Y"))

# check input
if not os.path.isfile(prm.POSCAR):
    print("\n** ERROR: Initial position file %s was not found." % prm.POSCAR)
    sys.exit(0)

if prm.get_super == True:
    # initialize Tmat
    Tmat = np.eye(3)
    # construct Tmat
    if len(prm.T_mat.split()) == 9:
        if all(check_int(s) for s in prm.T_mat.split()):
            #print(np.array(prm.T_mat.split(),dtype=int).reshape([3,3]))
            Tmat = np.array(prm.T_mat.split(),dtype=int).reshape([3,3])
            #Tmat = np.reshape(list(map(int, prm.T_mat.split()),[3,3]))
        else:
            print("\n** ERROR: Transformation matrix elements should be int.")
            sys.exit(0)
    elif len(list(map(int, prm.T_mat.split()))) == 3:
        if all(check_int(s) for s in prm.T_mat.split()):
            np.fill_diagonal(Tmat,list(map(int, prm.T_mat.split())))
        else:
            print("\n** ERROR: Transformation matrix elements should be int.")
            sys.exit(0)
    else:
        print("\n** ERROR: dimension of the transformation matrix wrong.")
        sys.exit(0)

if prm.getdist ==True:
    if len(prm.atom_list.split()) == 2:
        if all(check_int(s) for s in prm.atom_list.split()):
            atom_list = list(map(int, prm.atom_list.split()))
        else:
            print("\n** ERROR: Indexes of atoms should be int.")
            sys.exit(0)
    else:
        print("\n** ERROR: number indexes of atoms should be 2.")
        sys.exit(0)


# Read information from command line
# First specify location of POSCAR
i_POSCAR=prm.POSCAR.lstrip()

if prm.prt == True:
    print("\nPosition file name: %s " % i_POSCAR)
    # print("Symmetry tolerance: %s" % prm.symmetry)

# get cell informations
#----------------------------
initial_pos = read(i_POSCAR)#, format='vasp')
lattice     = initial_pos.get_cell()
positions   = initial_pos.get_scaled_positions()
numbers     = initial_pos.get_atomic_numbers()
# Magnetic moments are not considered in get_spacegroup method
#magmoms = [np.ones(initial_pos.get_number_of_atoms())]
#initial_cell = (lattice, positions, numbers)

if prm.prt == True:
    print('\n===========================================')
    print('\nInitial Structure')
    print("\nLattice Matrix  : (in Angstrom) ")
    print(np.array_str(lattice.array, precision=8))
    print("\nAtomic Positions: (in direct coordinate) ")
    print(positions)
    print("\nAtomic numbers  : (for each atom) ")
    print(numbers)
    print('\n===========================================')
# visualize
if prm.visualize == True:
    view(initial_pos)

# get atomic species
# fin = [initial_pos.get_chemical_symbols()[0]]
# for i in range(len(initial_pos.get_chemical_symbols())):
#     if i == len(initial_pos.get_chemical_symbols())-1: break
#     if initial_pos.get_chemical_symbols()[i] != initial_pos.get_chemical_symbols()[i+1]:
#         fin.append(initial_pos.get_chemical_symbols()[i+1])
atom_species=set(initial_pos.get_chemical_symbols()) # 2022-09-14: use set here.

# make supercell
#----------------------------
if prm.get_super == True:
    del initial_pos.constraints
    print(initial_pos,Tmat)
    super_cell = make_supercell(initial_pos,Tmat)

    # some post processing...
    #----------------------------
    indices = []
    symbls    = []
    for symbol in atom_species:
        ii = 0
        for atom in super_cell.get_chemical_symbols():
            if atom == symbol:
                indices.append(ii) # get sorted atomic index
                symbls.append(symbol) # get sorted atomic symbols
            ii += 1

    pos_new     = [super_cell.get_positions()[i] for i in indices]
    super_cell.set_positions(pos_new)
    super_cell.set_chemical_symbols(symbls)

    # get supercell informations
    #----------------------------
    lattice     = super_cell.get_cell()
    positions   = super_cell.get_scaled_positions()
    numbers     = super_cell.get_atomic_numbers()

    if prm.prt == True:
        print('\n===========================================')
        print('\nSuper Structure')
        print("\nLattice Matrix  : (in Angstrom) ")
        print(np.array_str(lattice.array, precision=8))
        print("\nAtomic Positions: (in direct coordinate) ")
        print(positions)
        print("\nAtomic numbers  : (for each atom) ")
        print(numbers)
        print('\n===========================================')
    # visualize?
    if prm.visualize==True:
        view(super_cell)

    # Write super structure to file
    if prm.over == True:
        write("POSCAR", super_cell, format='vasp')
    else:
        write(i_POSCAR+"_"+"super"+".vasp",super_cell,format='vasp')


    if float(sys.version.split()[0][:3]) < 3.0:
        # adding atomic species
        #----------------------------
        # get species
        # fin = [super_cell.get_chemical_symbols()[0]]
        # for i in range(len( super_cell.get_chemical_symbols())):
        #     if i == len( super_cell.get_chemical_symbols())-1: break
        #     if  super_cell.get_chemical_symbols()[i] !=  super_cell.get_chemical_symbols()[i+1]:
        #         fin.append( super_cell.get_chemical_symbols()[i+1])
        # atom_species=set(super_cell.get_chemical_symbols())

        # open file and readlines
        if prm.over == True:
            f = open("POSCAR", "r")
        else:
            f = open(i_POSCAR+"_"+"super"+".vasp", "r")
        contents = f.readlines()
        f.close()
        # format string
        contents.insert(5, ' '.join(atom_species)+'\n')
        # write string to file
        if prm.over == True:
            f = open("POSCAR", "r")
        else:
            f = open(i_POSCAR+"_"+"super"+".vasp", "w")
        contents = "".join(contents)
        f.write(contents)
        f.close()

    # output
    if prm.prt == True:
        if prm.over == False:
            print("\nOutput file name: %s " % str(i_POSCAR+"_"+"super"+".vasp"))


# get distance
#----------------------------
if prm.getdist == True:
    # here we use natural: atomic number starting from 1.
    dist = initial_pos.get_distance(atom_list[0]-1,atom_list[1]-1)
    if prm.prt == True:
        sys.stdout.write("\033[1;31m" ) # set color red
        print("\nDistance between atom %d and %d: %f" %(atom_list[0], atom_list[1], dist))
        sys.stdout.write("\033[0;0m") # reset color
    else:
        print(dist)


if prm.prt == True:
    # Post process
    #----------------------------
    endtime = time.time()
    runtime = endtime-starttime
    print("\nEnd of calculation.")
    print("Program was running for %.2f seconds." % runtime)
