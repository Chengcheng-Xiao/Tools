#!/bin/bash
# eliminate forces of fixed atoms
fix=$1
if test -z $1; then
fix=0
fi
# get force form OUTCAR
awk -v fix="$fix" '/POSITION/,/drift/{
if($1~/^[0-9.]+$/&&$3>=fix) print $1,$2,$3,sqrt($4*$4+$5*$5+$6*$6i);
else if($1=="total") print $0
}' OUTCAR >temp.f
awk '{
if($1=="total") {print ++i,a;a=0}
else {if(a<$4) a=$4}
}' temp.f >force.conv
#sed -i '1,9d' force.conv
#rm temp.f
tail -100 force.conv >temp.f
#get dist from XDATCAR
touch p1.conv p2.conv;
touch dist.conv
Num=`awk 'NR==7{for(i=1;i<=NF;i++) a=$i+a;}END{print a}' XDATCAR`
Lnum=`wc XDATCAR|awk '{print $1}'`
((n=(Lnum-7)/(Num+1)))
head -7 XDATCAR >p1.conv
awk -v num="$Num" 'NR==9,NR==(num+1)+7{print $0}' XDATCAR >>p1.conv
for((i=1;i<n;i++))
do
head -7 XDATCAR >p2.conv
((n1=9+(Num+1)*i))
((n2=(i+1)*(Num+1)+7))
sed -n ''$n1','$n2'p' XDATCAR >>p2.conv
echo -e $i"\t"`dist.pl p1.conv p2.conv ` >>dist.conv
done
#plot
gnuplot <<EOF
set term dumb
set xlabel 'Ion steps'
set ylabel 'Dist (Angstrom)'
plot 'dist.conv' w l t " Dist "
set xlabel 'Ion steps'
set ylabel 'Force (eV/Angstrom)'
plot 'temp.f' w l t " Force "
EOF
rm force.conv dist.conv p1.conv p2.conv temp.f