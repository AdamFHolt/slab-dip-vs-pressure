#!/bin/bash
# Progress of the running extract_properties.py suites.
#
# The per-model logs stay empty until each process exits (python buffers stdout
# when redirected), so progress is inferred from bytes read: each timestep
# reads exactly one full.N.csv, so rchar / csv_size is the timestep count.
cd /home/holt/Projects/ASPECT/subd_2D/compositional/analysis

declare -A MAXT=(
[2D_compositional_subd_lower-res_new_50plates]=25
[2D_compositional_subd_lower-res_new_FixedSP_50plates]=24
[2D_compositional_subd_lower-res_new_FixedOP_50plates]=25
[2D_compositional_subd_lower-res_new_250plates]=28
[2D_compositional_subd_lower-res_new_FixedSP_250plates]=30
[2D_compositional_subd_lower-res_new_FixedOP_250plates]=28
[2D_compositional_subd_lower-res_new2]=31
[2D_compositional_subd_FixedSP_lower-res_new2]=30
[2D_compositional_subd_FixedOP_lower-res_new]=28
[2D_compositional_subd_lower-res_new_1000plates]=30
[2D_compositional_subd_lower-res_new_FixedSP_1000plates2]=30
[2D_compositional_subd_lower-res_new_FixedOP_1000plates]=30
[2D_compositional_subd_lower-res_new_375plates]=30
[2D_compositional_subd_lower-res_new_FixedSP_375plates]=27
[2D_compositional_subd_lower-res_new_FixedOP_375plates]=32
)

printf '%-46s %-13s %7s %9s\n' MODEL OUTDIR ELAPSED PROGRESS
for pid in $(pgrep -f 'extract_properties\.py'); do
  [ -r "/proc/$pid/io" ] || continue
  # args follow the script name: <model> <maxt> <depth> <dz> <ds> <dz> <outdir> ...
  mapfile -t a < <(tr '\0' '\n' < "/proc/$pid/cmdline")
  i=0; for k in "${!a[@]}"; do [[ ${a[$k]} == *extract_properties.py ]] && i=$k; done
  model=${a[$((i+1))]}
  outdir=${a[$((i+7))]}
  [ -n "$model" ] || continue
  csv="csv_outputs/$model/full.20.csv"
  [ -f "$csv" ] || continue
  sz=$(stat -c %s "$csv")
  rc=$(awk '/^rchar/{print $2}' "/proc/$pid/io")
  done_t=$(( rc / sz ))
  tot=$(( ${MAXT[$model]:-30} - 8 ))
  el=$(ps -o etime= -p "$pid" | tr -d ' ')
  printf '%-46s %-13s %7s %4d/%-4d\n' \
      "${model#2D_compositional_subd_}" "$outdir" "$el" "$done_t" "$tot"
done | sort

echo
for d in TESTD_legacy TESTD_arc150; do
  n=$(ls text_files/$d/*.txt 2>/dev/null | grep -vc COLUMNS)
  echo "text_files/$d: ${n:-0}/15 models written"
done
