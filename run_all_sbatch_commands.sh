#!/bin/bash
set +x #echo on

for GCMfolder in GJ*; do
 cd ${GCMfolder}
 sbatch run_sbatch_greatlatkes
 cd ..
done

