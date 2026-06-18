#!/usr/bin/env bash

set -euo pipefail

for threads in 1 2 4 8 16; do
	export OMP_NUM_THREADS="${threads}"
	mkdir -p "Planet_Run_CK_par_${threads}"
	cp igcm3_nopg_par "Planet_Run_CK_par_${threads}/igcm3_nopg_par"
	cp fort.7corrk "Planet_Run_CK_par_${threads}/fort.7"
	cd "Planet_Run_CK_par_${threads}"
	ulimit -s unlimited
	echo "Running ./igcm3_nopg_par with OMP_NUM_THREADS=${threads}"
	start_time=$(date +%s.%N)
	OMP_NUM_THREADS="${threads}" ./igcm3_nopg_par >> runtime_par.txt 2>&1
	end_time=$(date +%s.%N)
	elapsed=$(echo "$end_time - $start_time" | bc)
	echo "Threads=${threads}: ${elapsed}s"
done
