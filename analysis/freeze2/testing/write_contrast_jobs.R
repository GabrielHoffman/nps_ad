# Gabriel Hoffman
# Write LSF job files to run dreamlet on contrasts

# Load CONTRASTS and metadata
load("/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/analysis/contrasts_for_dreamlet.Rdata")

write_job = function( variable_type, ctst_key, dataset, method){

	suffix = gsub(paste0(dataset, "_", dataset, "_"), paste0(dataset, "_"), paste0(toupper(dataset), "_", ctst_key, "_", SampleLevel, "_", AnnoLevel))

	outpath = paste0("/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/analysis/results/", suffix)

	outfile = paste0(outpath, "/", suffix, "_", method)

	if( ! dir.exists(outpath) ) dir.create(outpath)

	SRC = paste0("/sc/arion/projects/CommonMind/hoffman/NPS-AD/work/nps_ad/analysis/freeze2/testing/run_", method, "_contrasts.Rmd")

	cmd = paste0("rmarkdown::render('", SRC, "', 
output_file = '", outfile, "', 
params = list(dataset = '", dataset, "',
variable_type = '", variable_type, "', 
ctst_key = '", ctst_key, "'))")

	logs = paste0("/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/analysis/logs/", suffix, "_", method, "_")

	txt = paste0("#!/bin/bash
#BSUB -J ", suffix, "_", method, "
#BSUB -P acc_CommonMind
#BSUB -q premium
#BSUB -n 6
#BSUB -R rusage[mem=12000]
#BSUB -R span[hosts=1]
#BSUB -W 48:00
#BSUB -o ", logs, "%J.stdout
#BSUB -eo ", logs, "%J.stderr
#BSUB -L /bin/bash

ml hdf5/1.12.1 libpng/16 R/4.2.0 pandoc/2.6

R_LIBS_USER=/hpc/users/hoffmg01/.Rlib/R_420
R_LIBS=${R_LIBS_USER}:/hpc/users/hoffmg01/.Rlib/R_420:/hpc/packages/minerva-centos7/rpackages/4.2.0/site-library:/hpc/packages/minerva-centos7/rpackages/bioconductor/3.15

Rscript -e \"", cmd, "\"

")

	jobFile = paste0("/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/analysis/jobs/", suffix, "_", method, ".lsf") 

	write(txt, jobFile)
}


for( method in c("dreamlet", "crumblr") ){
	for( variable_type in names(CONTRASTS)){
		for( ctst_key in names(CONTRASTS[[variable_type]]) ){

			dataset = gsub("([a-zA-Z]+).*", "\\1", ctst_key)

			if( dataset %in% c("MSSM", "RUSH", "HBCC", "Aging") ){
				write_job( variable_type, ctst_key, dataset, method)
			}else{
				for(dataset in c("MSSM", "RUSH", "HBCC")){
					write_job( variable_type, ctst_key, dataset, method)
				}
			}
		}
	}
}

