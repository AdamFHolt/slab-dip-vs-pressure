#!/bin/bash
f=${1}

xmin=0;xmax=5800000;
ymin=0;ymax=1450000;

reg=-R$xmin/$xmax/$ymin/$ymax
int=-I2000
proj=-Jx0.000003

rm tmp.ps tmp1.grd tmp2.grd tmp3.grd slab.c.cpt input1.txt input2.txt input4.txt 2> /dev/null
gmt makecpt -Cpolar -T-0/1/0.001 -D > slab.c.cpt

echo 'converting text to grd files...'
gawk '{if(NR>2)print($1,$2,$3)}' outputs/"$f".txt > input1.txt
gawk '{if(NR>2)print($1,$2,$4)}' outputs/"$f".txt > input2.txt
gawk '{if(NR>2)print($1,$2,$5)}' outputs/"$f".txt > input3.txt

gmt surface $int $reg input1.txt -Gtmp1.grd
gmt surface $int $reg input2.txt -Gtmp2.grd
gmt surface $int $reg input3.txt -Gtmp3.grd

echo 'grd files created!'

gmt grdimage -Cslab.c.cpt tmp1.grd -Bn -K -Y+15 $reg $proj > tmp.ps
gmt grdimage -Cslab.c.cpt tmp2.grd -Bn -K -O -Y-5 $reg $proj >> tmp.ps
gmt grdimage -Cslab.c.cpt tmp3.grd -Bn -K -O -Y-5 $reg $proj >> tmp.ps
gmt psscale -D1/-1/2/.2h -Cslab.c.cpt -O -K -B0.5 >> tmp.ps

eps2eps tmp.ps tmp2.ps
mv tmp2.ps tmp.ps
convert -rotate 90 -flatten -trim -density 300 tmp.ps test_plots/$f.png

rm tmp.ps tmp1.grd tmp2.grd tmp3.grd slab.c.cpt input1.txt input2.txt input3.txt 2> /dev/null
