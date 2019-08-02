#!/bin/sh
#edit by lipai@USTC
#plot energies of structure in optimization
awk '/E0/{if ( i<=5 ) i++;else print $0 }' OSZICAR >temp.e
gnuplot <<EOF
set term dumb
set title 'Energy of each ion steps'
set xlabel 'Ion steps'
set ylabel 'Energy(eV)'
plot 'temp.e' u 1:5 w l
EOF
rm temp.e
