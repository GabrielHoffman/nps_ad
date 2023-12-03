#!/bin/bash 
#BSUB -J $1
#BSUB -P acc_CommonMind
#BSUB -q premium
#BSUB -n 12
#BSUB -R span[hosts=1]
#BSUB -R rusage[mem=50000]
#BSUB -W 96:00 
#BSUB -o /sc/arion/projects/CommonMind/hoffman/NPS-AD/work/nps_ad/analysis/freeze2/preprocess/logs/%J.stdout
#BSUB -eo /sc/arion/projects/CommonMind/hoffman/NPS-AD/work/nps_ad/analysis/freeze2/preprocess/logs/%J.stderr
#BSUB -L /bin/bash
#BSUB -cwd /sc/arion/projects/CommonMind/hoffman/NPS-AD/work/nps_ad/analysis/freeze2/preprocess/

source /hpc/users/hoffmg01/.bash_profile
ml R/4.3.0

echo $1
./submit_preproc.R --cohort $1

