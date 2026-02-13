#!/bin/bash

m="2D_compositional_subd_lower-res_new_FixedSP_50plates"
maxt="28"

depth="300.0e3"

python3 extract_properties_and_plot2.py $m $maxt $depth 10.e3 10.0e3 1.0e3

depth="400.0e3"

python3 extract_properties_and_plot2.py $m $maxt $depth 10.e3 10.0e3 1.0e3






























