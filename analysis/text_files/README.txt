What each directory here holds.

  extracted_archive/    The extraction the paper's figures rest on. 21 columns.
                        Written by an older extract_properties.py, so the
                        current code does not reproduce its stress columns
                        exactly (geometry it does, to machine precision).
                        Treat it as a frozen archive.

  extracted_Lscales/    Default output of the current extract_properties.py.
                        Adds the velocity and length-scale diagnostics used by
                        plot_Leff_distribution.supp.py (v_s, dv_s/ds, L_v, L_K,
                        the dQ/ds scaling, slab viscosity).

  extracted_coverage/   Re-extraction at 250, 300 and 350 km carrying the
                        midplane coverage columns 30-38 (see COLUMNS.txt
                        beside it). Built by many_property-extractions.recheck.sh.

  withT/                extracted_archive with column 17 upgraded from the
                        deviatoric normal-stress term H*tau_n_bar to the full
                        curvature contribution K*H*(sigma_n_bar - P_ext).
                        Built by make_withT.py. Every force figure reads this.

  withT_covered/        withT with the timesteps that fail the coverage test
                        blanked to NaN. Built by make_withT_covered.py, read
                        by the *.covered.py figure variants.

  curvature_pressure_term/  The measured pressure part T, per model and depth,
                        plus the SUMMARY*.txt files that back the numbers
                        quoted in the manuscript.

  Pref/                 Far-field reference pressure from extract_pref.py,
                        read by plot_Psp-DP-ratio.py.

Only the SUMMARY*.txt and BALANCE_COMPARISON.txt files are tracked in git.
Everything else regenerates from the scripts in the parent directory.
