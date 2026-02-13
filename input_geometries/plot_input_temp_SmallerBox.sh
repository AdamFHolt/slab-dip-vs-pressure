#!/bin/bash
f=${1}

xmin=0;xmax=5800000;
ymin=0;ymax=1450000;

reg=-R$xmin/$xmax/$ymin/$ymax
int=-I2000
proj=-Jx0.000003

rm tmp.ps tmp.grd slab.c.cpt input.txt 2> /dev/null

gawk '{if(NR>2)print($1,$2,$3)}' outputs/"$f".txt > input.txt

gmt makecpt -Cpolar -T-200/1600/0.001 -D > slab.c.cpt

echo 'converting text to grd file...'
gmt surface $int $reg input.txt -Gtmp.grd
echo 'grd file created!'

gmt grdimage -Cslab.c.cpt tmp.grd -Ba1000000:"x":/a1000000:"y":WeSn -K $reg $proj > tmp.ps
gmt psscale -D1/-1/5/.2h -Cslab.c.cpt -O -K -B500 >> tmp.ps

eps2eps tmp.ps tmp2.ps
mv tmp2.ps tmp.ps
convert -rotate 90 -flatten -trim -density 300 tmp.ps test_plots/$f.png
echo written to slabslices/$f.png

rm tmp.ps tmp.grd slab.c.cpt input.txt 2> /dev/null
