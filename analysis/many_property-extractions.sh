#!/bin/bash

depth="300.0e3"
script="extract_properties.py"
MAX_JOBS="${MAX_JOBS:-8}"

run_job () {
    local model="$1"
    local maxt="$2"
    echo "[START] ${model} (maxt=${maxt})"
    python3 "$script" "$model" "$maxt" "$depth" 10.0e3 10.0e3 1.0e3
    local rc=$?
    if [ $rc -ne 0 ]; then
        echo "[FAIL]  ${model} (exit=${rc})"
    else
        echo "[DONE]  ${model}"
    fi
    return $rc
}

queue_job () {
    run_job "$1" "$2" &
    while [ "$(jobs -rp | wc -l)" -ge "$MAX_JOBS" ]; do
        wait -n
    done
}

#--

m="2D_compositional_subd_lower-res_new_50plates"
maxt=25
queue_job "$m" "$maxt"

m="2D_compositional_subd_lower-res_new_FixedSP_50plates"
maxt=24
queue_job "$m" "$maxt"

m="2D_compositional_subd_lower-res_new_FixedOP_50plates"
maxt=25
queue_job "$m" "$maxt"

#--

m="2D_compositional_subd_lower-res_new_250plates"
maxt=28
queue_job "$m" "$maxt"
	
m="2D_compositional_subd_lower-res_new_FixedSP_250plates"
maxt=30
queue_job "$m" "$maxt"

m="2D_compositional_subd_lower-res_new_FixedOP_250plates"
maxt=28
queue_job "$m" "$maxt"

#--

m="2D_compositional_subd_lower-res_new2"      
maxt=31
queue_job "$m" "$maxt"

m="2D_compositional_subd_FixedSP_lower-res_new2"
maxt=30
queue_job "$m" "$maxt"

m="2D_compositional_subd_FixedOP_lower-res_new"
maxt=28
queue_job "$m" "$maxt"

#--

m="2D_compositional_subd_lower-res_new_1000plates"
maxt=30
queue_job "$m" "$maxt"

m="2D_compositional_subd_lower-res_new_FixedSP_1000plates2"
maxt=30
queue_job "$m" "$maxt"

m="2D_compositional_subd_lower-res_new_FixedOP_1000plates"
maxt=30
queue_job "$m" "$maxt"

#--

m="2D_compositional_subd_lower-res_new_375plates"
maxt=30
queue_job "$m" "$maxt"

m="2D_compositional_subd_lower-res_new_FixedSP_375plates"
maxt=27
queue_job "$m" "$maxt"

m="2D_compositional_subd_lower-res_new_FixedOP_375plates"
maxt=32
queue_job "$m" "$maxt"

wait
echo "All extraction jobs finished."
