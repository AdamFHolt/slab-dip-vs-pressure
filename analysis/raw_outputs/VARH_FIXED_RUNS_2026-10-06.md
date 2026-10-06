# Variable-H runs, second set (2026-10-06): fixed-SP, fixed-OP, and one B-matched run

13 runs. All prm files are the corresponding free-plate H prm (new_{250,500,1000}plates_H{70,130}km,
stampede3 paths, End time 90 Myr, checkpoint every 25 steps) with only the output directory and the
composition file changed; the drho71 run also changes the lithosphere density.

## Files to copy to the cluster

prm and .batch files (this directory; one .batch per run, 2 nodes / 96 tasks / skx / 48 h, same recipe as 500-130km) -> /work2/04714/adamholt/stampede3/aspect_work/subd_models_new/compositional/
  2D_compositional_subd_lower-res_new_FixedSP_{250,500,1000}plates_H{70,130}km.prm
  2D_compositional_subd_lower-res_new_FixedOP_{250,500,1000}plates_H{70,130}km.prm
  2D_compositional_subd_lower-res_new_500plates_H70km_drho71.prm

composition files (input_geometries/outputs/, gzipped, ~22 MB each) -> .../compositional/text_files/
  comp_lith-and-crust_FixedSP_H{70,130}km.txt.gz
  comp_lith-and-crust_FixedOP_tapered_H{70,130}km.txt.gz
  (gunzip on the cluster; the drho71 run reuses comp_lith-and-crust_tapered_H70km.txt, already there)

Generators: input_geometries/make_comp_lith-and-crust_FixedSP_H.py and
make_comp_lith-and-crust_FixedOP_tapered_H.py <H_km>; both reproduce the original H = 100 files byte for byte.

## Suggested job names and wall times

The 500 and 1000 free-plate runs at H = 130 used 32 h / 48 h and the 500 one ran out at 89.6 Myr;
the analysis window needs ~60 Myr (slab tip past ~1250 km), so 90 Myr is not required.

  SP-250-70   SP-500-70   SP-1000-70      48 h each
  SP-250-130  SP-500-130  SP-1000-130     48 h
  OP-250-70   OP-500-70   OP-1000-70      48 h
  OP-250-130  OP-500-130  OP-1000-130     48 h   (500 and 1000 expected to overturn; kept deliberately)
  drho71-500-70                           48 h

## Notes for the analysis

- The drho71 run holds B = drho*g*H = 71.43*9.81*70 km = 49.05 MPa (same as H = 100). extract_properties.py
  has drho = 50 hard-coded (DP_anal); that run needs drho = 71.43 passed in before extraction.
- nominal_plate_thickness() parses _H(\d+)km, so all 13 names work with the varH scripts as they are.
- Expected classification: fixed-SP all normal (never overturned in the suite); fixed-OP 500/1000 at
  H = 130 overturned; fixed-OP at H = 70 probably normal.

## Submit

  cd /work2/04714/adamholt/stampede3/aspect_work/subd_models_new/compositional
  for f in *_Fixed??_*_H*km.batch *_drho71.batch; do sbatch $f; done
