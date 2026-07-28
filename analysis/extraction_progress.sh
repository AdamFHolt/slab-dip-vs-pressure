#!/bin/bash
# Progress of any running extract_properties.py suites.
#
# The per-model logs stay empty until each process exits (python buffers stdout
# when redirected), so progress is inferred from bytes read: each timestep
# reads exactly one full.N.csv, so rchar / csv_size is the timestep count.
# Completed models only appear in the tally at the bottom once they finish,
# because np.savetxt runs once at the end.
#
# Usage:  bash extraction_progress.sh
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

tmp=$(mktemp); nrun=0
for pid in $(pgrep -f 'extract_properties\.py'); do
  [ -r "/proc/$pid/io" ] || continue
  # args follow the script name: <model> <maxt> <depth> <dz> <ds> <dz> <outdir> ...
  mapfile -t a < <(tr '\0' '\n' < "/proc/$pid/cmdline")
  i=0; for k in "${!a[@]}"; do [[ ${a[$k]} == *extract_properties.py ]] && i=$k; done
  model=${a[$((i+1))]}; depth=${a[$((i+3))]}; outdir=${a[$((i+7))]}
  [ -n "$model" ] || continue
  csv="csv_outputs/$model/full.20.csv"
  [ -f "$csv" ] || continue
  sz=$(stat -c %s "$csv")
  rc=$(awk '/^rchar/{print $2}' "/proc/$pid/io")
  printf '%-46s %-13s %6s %9s %4d/%-4d\n' \
      "${model#2D_compositional_subd_}" "$outdir" \
      "$(awk -v d="$depth" 'BEGIN{printf "%.0fkm", d/1000}')" \
      "$(ps -o etime= -p "$pid" | tr -d ' ')" \
      "$(( rc / sz ))" "$(( ${MAXT[$model]:-30} - 8 ))" >> "$tmp"
  nrun=$((nrun+1))
done

if [ "$nrun" -gt 0 ]; then
  printf '%-46s %-13s %6s %9s %9s\n' MODEL OUTDIR DEPTH ELAPSED PROGRESS
  sort -k3,3 -k1,1 "$tmp"
else
  echo "no extract_properties.py processes running"
fi
rm -f "$tmp"

echo
echo "completed models, per output directory and depth:"
for d in text_files/TESTD_*; do
  [ -d "$d" ] || continue
  for z in $(ls "$d"/*.txt 2>/dev/null | grep -o '\.z[0-9.]*\.shear' \
             | sed 's/^\.//;s/\.shear$//' | sort -u); do
    n=$(ls "$d"/*".$z."*.txt 2>/dev/null | wc -l)
    w=$(cat "$d"/logs/*".$z.log" 2>/dev/null | grep -c '^WARNING')
    printf '  %-16s %-8s %2d/15 models   %2d timesteps flagged as not covered\n' \
        "${d#text_files/}" "$z" "$n" "$w"
  done
done
