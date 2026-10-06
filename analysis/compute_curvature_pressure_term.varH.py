"""
Run compute_curvature_pressure_term.py on the variable-H comparison set without
touching its default inputs or outputs.

The six H = 70/130 km runs are read from text_files/extracted_varH/ and their
three H = 100 counterparts from text_files/extracted_coverage/, both written by
the current extract_properties.py with the same arguments.  (The paper's
curvature_pressure_term/ files were built from extracted_archive/, an older
extraction, so they are not reused here.)  Output goes to
text_files/curvature_pressure_term_varH/.

Usage:  python3 compute_curvature_pressure_term.varH.py <model> [depth_m]
        bash many_curvature-pressure-term.varH.sh [depth_m]
"""
import sys

import compute_curvature_pressure_term as cpt

H100 = ('2D_compositional_subd_lower-res_new_250plates',
        '2D_compositional_subd_lower-res_new2',
        '2D_compositional_subd_lower-res_new_1000plates')

src = 'extracted_coverage' if sys.argv[1] in H100 else 'extracted_varH'
cpt.TXT = cpt.TXT.replace('extracted_archive', src)
cpt.OUTDIR = cpt.ANALYSIS + 'text_files/curvature_pressure_term_varH'
cpt.main()
