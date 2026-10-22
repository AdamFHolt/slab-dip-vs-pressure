"""
Run compute_curvature_pressure_term.py on the early outputs (t = 5-8) of
text_files/extracted_early/ (many_property-extractions.early.sh), for the
tip-depth test of reviewer 2, comment 3.  Output goes to
text_files/curvature_pressure_term_early/; the default inputs and outputs are
untouched.

Usage:  python3 compute_curvature_pressure_term.early.py <model> [depth_m]
"""
import compute_curvature_pressure_term as cpt

cpt.TXT = cpt.TXT.replace('extracted_archive', 'extracted_early')
cpt.OUTDIR = cpt.ANALYSIS + 'text_files/curvature_pressure_term_early'
cpt.first_time = 5
cpt.start_time = 5
cpt.main()
