#!/bin/bash

#$ -cwd

#$ -j y

#$ -S /bin/bash

#$ -pe smp 2

#$ -l h_vmem=8g

#$ -l h_rt=24:00:00


##export OMP_NUM_THREADS=1
/share/apps/matlab_current/bin/matlab -r "try; Crunch_Run_Timing($M, $pqi_mode, $var_mode); catch me; display(me); end; quit"

## call with: qsub -v M=1000,pqi_mode=1,var_mode=1 run_script_timings.sh