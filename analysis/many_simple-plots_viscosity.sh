#!/bin/bash

script="plot_onestep_simple-viscosity.zoomed.py"


m1="2D_compositional_subd_lower-res_new_50plates"
t1=16
x1=3250
m2="2D_compositional_subd_lower-res_new_FixedSP_50plates"
t2=16
x2=3200
m3="2D_compositional_subd_lower-res_new_FixedOP_50plates"
t3=16
x3=3250
m4="2D_compositional_subd_lower-res_new_250plates"
t4=16
x4=3250
m5="2D_compositional_subd_lower-res_new_FixedSP_250plates"
t5=16
x5=3000
m6="2D_compositional_subd_lower-res_new_FixedOP_250plates"
t6=16
x6=3350
m7="2D_compositional_subd_lower-res_new2"         # 500
t7=17
x7=3350
m8="2D_compositional_subd_FixedSP_lower-res_new2" # 500 
t8=17
x8=2900
m9="2D_compositional_subd_FixedOP_lower-res_new"  # 500
t9=17
x9=3350
m10="2D_compositional_subd_lower-res_new_1000plates"
t10=18
x10=3250
m11="2D_compositional_subd_lower-res_new_FixedSP_1000plates2"
t11=18
x11=2900
m12="2D_compositional_subd_lower-res_new_FixedOP_1000plates"
t12=18
x12=3350
m13="2D_compositional_subd_lower-res_new_375plates"
t13=17
x13=3250
m14="2D_compositional_subd_lower-res_new_FixedSP_375plates"
t14=17
x14=3000
m15="2D_compositional_subd_lower-res_new_FixedOP_375plates"
t15=17
x15=3350


python $script $m1 $t1 $x1 &
python $script $m2 $t2 $x2 &
python $script $m3 $t3 $x3 &
python $script $m4 $t4 $x4 &
python $script $m5 $t5 $x5 &
python $script $m6 $t6 $x6 &
python $script $m7 $t7 $x7 &
python $script $m8 $t8 $x8 &
python $script $m9 $t9 $x9 &
python $script $m10 $t10 $x10 &
python $script $m11 $t11 $x11 &
python $script $m12 $t12 $x12 &
python $script $m13 $t13 $x13 &
python $script $m14 $t14 $x14 &
python $script $m15 $t15 $x15 &

wait
echo "All plots done"
