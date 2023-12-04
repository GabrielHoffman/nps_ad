#!/bin/bash 
#BSUB -J $1
#BSUB -P acc_CommonMind
#BSUB -q premium
#BSUB -n 12
#BSUB -R span[hosts=1]
#BSUB -R rusage[mem=50000]
#BSUB -W 96:00 
#BSUB -o /sc/arion/projects/CommonMind/hoffman/NPS-AD/work/nps_ad/analysis/freeze2/preprocess/logs/$1_%J.stdout
#BSUB -eo /sc/arion/projects/CommonMind/hoffman/NPS-AD/work/nps_ad/analysis/freeze2/preprocess/logs/$1_%J.stderr
#BSUB -L /bin/bash
#BSUB -cwd /sc/arion/projects/CommonMind/hoffman/NPS-AD/work/nps_ad/analysis/freeze2/preprocess/

source /hpc/users/hoffmg01/.bash_profile
ml R/4.3.0 pandoc

export R_LIBS=/hpc/users/hoffmg01/.Rlib/R_430:/hpc/packages/minerva-centos7/rpackages/4.3.0/site-library:/hpc/packages/minerva-centos7/rpackages/bioconductor/3.17

echo $1
./submit_preproc.R --cohort $1

