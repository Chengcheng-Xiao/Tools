#!/usr/bin/env python
"""
A script to replace some naming convention in vasprun.xml


v5.4.4    v5.4.1

x2-y2 -->   dx2
fy3x2 -->    f1
 fxyz -->    f2
 fyz2 -->    f3
  fz3 -->    f4
 fxz2 -->    f5
 fzx2 -->    f6
  fx3 -->    f7

usage ./p4vmod.py vasprun.xml
"""
import sys
import os

# Find out how many arguments were on the command line,
nsubtract = len(sys.argv)-1
if nsubtract > 1:
    print "\n** ERROR: Must specify only one file."
    print "eg. ./chgdiff.py vasprun.xml"
    sys.exit(0)

elif nsubtract == 1 and not os.path.isfile(sys.argv[1]):
    print "\n** ERROR: Input file %s was not found." % sys.argv[1]
    sys.exit(0)

elif nsubtract == 1 and os.path.isfile(sys.argv[1]):
    filename = sys.argv[1].lstrip()

elif nsubtract == 0:
    print "\n** Warning: Defult to vasprun.xml."
    filename = "vasprun.xml"


# open file
f = open(filename, "r")
contents = f.readlines()
f.close()

# concatenat
contents = "".join(contents)

# replace
contents = contents.replace("x2-y2","  dx2")
contents = contents.replace("fy3x2","   f1")
contents = contents.replace(" fxyz","   f2")
contents = contents.replace(" fyz2","   f3")
contents = contents.replace("  fz3","   f4")
contents = contents.replace(" fxz2","   f5")
contents = contents.replace(" fzx2","   f6")
contents = contents.replace("  fx3","   f7")

# output
f = open(filename+".541.xml", "w")
f.write(contents)
f.close()
