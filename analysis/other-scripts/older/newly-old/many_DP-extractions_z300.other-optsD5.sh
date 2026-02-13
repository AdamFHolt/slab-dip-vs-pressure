#!/bin/bash

depth="300.0e3"
script="extract_properties_and_plot.simplified.py"



m="2D_compositional_subd_lower-res_new_375plates"
maxt=30
python3 $script $m $maxt $depth 5.0e3 10.e3 1.e3

m="2D_compositional_subd_lower-res_new_FixedSP_375plates"
maxt=27
python3 $script $m $maxt $depth 5.0e3 10.e3 1.e3

m="2D_compositional_subd_lower-res_new_FixedOP_375plates"
maxt=32
python3 $script $m $maxt $depth 5.0e3 10.e3 1.e3


