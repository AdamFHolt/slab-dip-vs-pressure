#!/bin/bash

depth="300.0e3"
script="extract_properties_and_plot.simplified.py"


#--

m="2D_compositional_subd_lower-res_new_1000plates"
maxt=30
python3 $script $m $maxt $depth 5.0e3 10.e3 1.e3

m="2D_compositional_subd_lower-res_new_FixedSP_1000plates2"
maxt=30
python3 $script $m $maxt $depth 5.0e3 10.e3 1.e3

m="2D_compositional_subd_lower-res_new_FixedOP_1000plates"
maxt=30
python3 $script $m $maxt $depth 5.0e3 10.e3 1.e3

