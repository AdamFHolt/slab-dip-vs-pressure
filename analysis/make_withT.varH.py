"""
Build text_files/withT_varH/ for the variable-H comparison set by running
make_withT.py with its inputs and output redirected; make_withT.py itself and
text_files/withT/ are untouched.

  H = 70/130 km runs:  extracted_varH/     + curvature_pressure_term_varH/
  H = 100 km runs:     extracted_coverage/ + curvature_pressure_term_varH/

Run many_property-extractions.varH.sh and many_curvature-pressure-term.varH.sh
(at 250, 300 and 350 km) first.

Usage:  python3 make_withT.varH.py [depth_m ...]  (default 250e3 300e3 350e3)
"""
import os

import make_withT as mw

TF = os.path.join(mw.HERE, 'text_files')
mw.DST = os.path.join(TF, 'withT_varH')
mw.TDIR = os.path.join(TF, 'curvature_pressure_term_varH')

SETS = [
    ('extracted_coverage', ['2D_compositional_subd_lower-res_new_250plates',
                            '2D_compositional_subd_lower-res_new2',
                            '2D_compositional_subd_lower-res_new_1000plates']),
    ('extracted_varH', ['2D_compositional_subd_lower-res_new_%dplates_H%dkm' % (e, h)
                        for e in (250, 500, 1000) for h in (70, 130)]),
]

for src, models in SETS:
    mw.SRC = os.path.join(TF, src)
    mw.MODELS = models
    mw.main()
