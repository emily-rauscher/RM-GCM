#!/bin/bash
set +x #echo on

for GCMfolder in HD*; do
 cd ${GCMfolder}
 sbatch run_sbatch_greatlatkes_restart
 cd ..
done

