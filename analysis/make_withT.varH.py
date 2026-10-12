"""
Build text_files/withT_varH/ for the variable-H comparison set by running
make_withT.py with its inputs and output redirected; make_withT.py itself and
text_files/withT/ are untouched.

  H = 70/130 km runs:  extracted_varH/     + curvature_pressure_term_varH/
  H = 100 km runs:     extracted_coverage/ + curvature_pressure_term_varH/

Free, fixed-SP and fixed-OP plates, plus the B-matched drho71 run.  Models with
no extraction file yet are skipped (the fixed-plate runs come in over time).

Run many_property-extractions.varH.sh and many_curvature-pressure-term.varH.sh
(at 250, 300 and 350 km) first.

Usage:  python3 make_withT.varH.py [depth_m ...]  (default 250e3 300e3 350e3)
"""
import os

import make_withT as mw

TF = os.path.join(mw.HERE, 'text_files')
mw.DST = os.path.join(TF, 'withT_varH')
mw.TDIR = os.path.join(TF, 'curvature_pressure_term_varH')

H100 = {'': ['2D_compositional_subd_lower-res_new_250plates',
             '2D_compositional_subd_lower-res_new2',
             '2D_compositional_subd_lower-res_new_1000plates'],
        'FixedSP': ['2D_compositional_subd_lower-res_new_FixedSP_250plates',
                    '2D_compositional_subd_FixedSP_lower-res_new2',
                    '2D_compositional_subd_lower-res_new_FixedSP_1000plates2'],
        'FixedOP': ['2D_compositional_subd_lower-res_new_FixedOP_250plates',
                    '2D_compositional_subd_FixedOP_lower-res_new',
                    '2D_compositional_subd_lower-res_new_FixedOP_1000plates']}
VARH = (['2D_compositional_subd_lower-res_new_%s%dplates_H%dkm' % (bc, e, h)
         for bc in ('', 'FixedSP_', 'FixedOP_') for e in (250, 500, 1000) for h in (70, 130)]
        + ['2D_compositional_subd_lower-res_new_500plates_H70km_drho71'])

SETS = [('extracted_coverage', sum(H100.values(), [])),
        ('extracted_varH', VARH)]

def have(src, m):
    return os.path.exists(os.path.join(TF, src, m + '.z300.0' + mw.SUFFIX))

for src, models in SETS:
    missing = [m for m in models if not have(src, m)]
    for m in missing:
        print('skipping %s (not extracted yet)' % m)
    mw.SRC = os.path.join(TF, src)
    mw.MODELS = [m for m in models if m not in missing]
    mw.main()
