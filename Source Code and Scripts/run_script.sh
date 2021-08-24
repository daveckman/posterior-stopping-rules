#!/bin/bash

#$ -cwd

#$ -j y

#$ -S /bin/bash

#$ -pe smp 8

#$ -l h_vmem=8g

#$ -l h_rt=24:00:00


##export OMP_NUM_THREADS=1
/share/apps/matlab_current/bin/matlab -r "try; Crunch_Run_Efficiency_Experiment_DoubleSplit($M, $Q, $pqi_mode, $rpi_mode); catch me; display(me); end; quit"

## call with: qsub -v M=100,Q=50,pqi_mode=1,rpi_mode=1 run_script.sh