#!/bin/bash

depth="250.0e3"
script="extract_properties_and_plot.simplified.py"

#--

m="2D_compositional_subd_lower-res_new_50plates"
maxt=27
python3 $script $m $maxt $depth 10.e3 10.0e3 1.0e3

m="2D_compositional_subd_lower-res_new_FixedSP_50plates"
maxt=28
python3 $script $m $maxt $depth 10.e3 10.0e3 1.0e3

m="2D_compositional_subd_lower-res_new_FixedOP_50plates"
maxt=27
python3 $script $m $maxt $depth 10.e3 10.0e3 1.0e3

#--

m="2D_compositional_subd_lower-res_new_250plates"
maxt=28
python3 $script $m $maxt $depth 10.e3 10.0e3 1.0e3
	
m="2D_compositional_subd_lower-res_new_FixedSP_250plates"
maxt=30
python3 $script $m $maxt $depth 10.e3 10.0e3 1.0e3

m="2D_compositional_subd_lower-res_new_FixedOP_250plates"
maxt=28
python3 $script $m $maxt $depth 10.e3 10.0e3 1.0e3

#--

m="2D_compositional_subd_lower-res_new2"      
maxt=31
python3 $script $m $maxt $depth 10.e3 10.0e3 1.0e3

m="2D_compositional_subd_FixedSP_lower-res_new2"
maxt=30
python3 $script $m $maxt $depth 10.e3 10.0e3 1.0e3

m="2D_compositional_subd_FixedOP_lower-res_new"
maxt=28
python3 $script $m $maxt $depth 10.e3 10.0e3 1.0e3

#--

m="2D_compositional_subd_lower-res_new_1000plates"
maxt=30
python3 $script $m $maxt $depth 10.e3 10.0e3 1.0e3

m="2D_compositional_subd_lower-res_new_FixedSP_1000plates2"
maxt=30
python3 $script $m $maxt $depth 10.e3 10.0e3 1.0e3

m="2D_compositional_subd_lower-res_new_FixedOP_1000plates"
maxt=30
python3 $script $m $maxt $depth 10.e3 10.0e3 1.0e3

#--

m="2D_compositional_subd_lower-res_new_375plates"
maxt=30
python3 $script $m $maxt $depth 10.e3 10.0e3 1.0e3

m="2D_compositional_subd_lower-res_new_FixedSP_375plates"
maxt=27
python3 $script $m $maxt $depth 10.e3 10.0e3 1.0e3

m="2D_compositional_subd_lower-res_new_FixedOP_375plates"
maxt=32
python3 $script $m $maxt $depth 10.e3 10.0e3 1.0e3


