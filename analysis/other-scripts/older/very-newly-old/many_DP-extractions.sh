#!/bin/bash

depth="300.0e3"
script="extract_properties.py"

#--

m="2D_compositional_subd_lower-res_new_50plates"
maxt=27
python3 $script $m $maxt $depth 10.0e3 10.0e3 1.0e3

m="2D_compositional_subd_lower-res_new_FixedSP_50plates"
maxt=28
python3 $script $m $maxt $depth 10.0e3 10.0e3 1.0e3

m="2D_compositional_subd_lower-res_new_FixedOP_50plates"
maxt=27
python3 $script $m $maxt $depth 10.0e3 10.0e3 1.0e3



