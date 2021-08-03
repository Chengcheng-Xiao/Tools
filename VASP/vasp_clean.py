#!/usr/bin/env python
"""
A script to clean current working dir.
"""
import os
import argparse

# Command line praser
#----------------------------
parser = argparse.ArgumentParser(description='A script to clean current working dir.')
parser.add_argument('-f','--force', action="store_true", default=False, dest='force',
                    help='Remove without warning? (Default=False)')
parser.add_argument('-a','--add', nargs='+', action="store", default=None, dest="additional",
                    help='Additional files to retain.')

prm = parser.parse_args()

# get a list of objects under current dir
dir_file = os.listdir('./')

# target files
files = ['INCAR', 'KPOINTS', 'POSCAR', 'POTCAR', 'vdw_kernel.bindat', '.pbs', 'wanier90']

if prm.additional != None:
    files += prm.additional

# exclude target files
for i in files:
    dir_file = [x for x in dir_file if not x.startswith(i)]
    dir_file = [x for x in dir_file if not x.endswith(i)]

# exclude directories
dir_file = [i for i in dir_file if not os.path.isdir(i)]

if prm.force == False:
    # Warning is very much needed!
    print 'Removing:\n\t','\n\t'.join(dir_file)
    val = raw_input("Are you sure? (y/n) ")
    if val == 'y':
        # remove files
        for file in dir_file:
            os.remove(file)
        print 'done.'
    elif val == 'n':
        print "exiting..."
    else:
        print('Input not accepted, exiting...')
        exit()
else:
    # remove files
    for file in dir_file:
        os.remove(file)

