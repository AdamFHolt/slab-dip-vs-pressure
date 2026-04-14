#!/bin/bash

# H, curvature, vc maps
python3 plot_H-curvature-vc_maps.py

# dip, age maps
python3 plot_dip-age_maps.py

# just-maps: ref, low visc, high visc
python3 plot_final_just-maps.py 4e22 3.28e-5 1333 8.044e-7 88e3
python3 plot_final_just-maps.py 2e22 3.28e-5 1333 8.044e-7 88e3
python3 plot_final_just-maps.py 8e22 3.28e-5 1333 8.044e-7 88e3
