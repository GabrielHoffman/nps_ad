#!/bin/bash
#BSUB -J $1
#BSUB -P acc_CommonMind
#BSUB -q premium
#BSUB -n 12
#BSUB -R span[hosts=1]
#BSUB -R usage[mem=50000]
#BSUB -W 96:00 
#BSUB -o ./${1}_%J.stdout
#BSUB -eo ./${1}_%J.stderr
#BSUB -L /bin/bash
#BSUB -cwd /sc/arion/projects/CommonMind/hoffman/NPS-AD/work/nps_ad/analysis/freeze2/preprocess/

ml purge
ml python git pandoc gcc/11.2.0
ml R/4.3.0 
export R_LIBS_USER=/hpc/users/hoffmg01/.Rlib/R_430/
export R_LIBS=\$R_LIBS_USER:\$R_LIBS

./submit_preproc.R --cohort $1

