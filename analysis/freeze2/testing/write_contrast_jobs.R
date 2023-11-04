# Gabriel Hoffman
# Write LSF job files to run dreamlet on contrasts

# Load CONTRASTS and metadata
# synapse: syn51114763
load("/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/analysis/contrasts_for_dreamlet.Rdata")

write_job = function( variable_type, ctst_key, dataset, method, SampleLevel, AnnoLevel){

	suffix = gsub(paste0(dataset, "_", dataset, "_"), paste0(dataset, "_"), paste0(toupper(dataset), "_", ctst_key, "_", SampleLevel, "_", AnnoLevel))

	suffix2 = gsub(paste0(dataset, "/", dataset, "/"), paste0(dataset, "/"), paste0(toupper(dataset), "/", ctst_key, "/", SampleLevel, "/", AnnoLevel))

	outpath = paste0("/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/analysis/results/", suffix2)

	outfile = paste0(outpath, "/", suffix, "_", method)

	if( ! dir.exists(outpath) ) dir.create(outpath, recursive=TRUE)

	SRC = paste0("/sc/arion/projects/CommonMind/hoffman/NPS-AD/work/nps_ad/analysis/freeze2/testing/run_", method, "_contrasts.Rmd")

	cmd = paste0("rmarkdown::render('", SRC, "', 
output_file = '", outfile, "',  
params = list(dataset = '", dataset, "',
variable_type = '", variable_type, "', 
SampleLevel = '", SampleLevel, "', 
AnnoLevel = '", AnnoLevel, "', 
ctst_key = '", ctst_key, "'))")

	logs = paste0("/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/analysis/logs/", suffix, "_", method, "_")

	txt = paste0("#!/bin/bash
#BSUB -J ", suffix, "_", method, "
#BSUB -P acc_CommonMind
#BSUB -q premium
#BSUB -n 10
#BSUB -R rusage[mem=16000]
#BSUB -R span[hosts=1]
#BSUB -W 48:00
#BSUB -o ", logs, "%J.stdout
#BSUB -eo ", logs, "%J.stderr
#BSUB -L /bin/bash

cd ",
outpath, 
"
ml hdf5/1.12.1 libpng/16 R/4.3.0 pandoc/2.6

R_LIBS_USER=/hpc/users/hoffmg01/.Rlib/R_430
R_LIBS=${R_LIBS_USER}:/hpc/packages/minerva-centos7/rpackages/4.3.0/site-library:/hpc/packages/minerva-centos7/rpackages/bioconductor/3.17

Rscript -e \"", cmd, "\"

")

	jobFile = paste0("/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/analysis/jobs/", suffix, "_", method, ".lsf") 

	write(txt, jobFile)
}

# create jobs
AnnoLevel = c("bulk", "class", "subclass", "subtype")
SampleLevel = c("Channel", "SubID")
method = c("dreamlet", "crumblr")
# dataset = c("MSSM", "RUSH", "HBCC", "AGING", "FULL")

df = expand.grid( 	SampleLevel = SampleLevel, 
					AnnoLevel = AnnoLevel, 
					variable_type = names(CONTRASTS),
					method = method)

# Remove bulk/dreamlet
df = df[!with(df, AnnoLevel=="bulk" & method == 'crumblr'),]

for(i in seq(nrow(df))){

	for(ctst_key in names(CONTRASTS[[df$variable_type[i]]]) ){

		dataset = gsub("([a-zA-Z]+).*", "\\1", ctst_key)

		if( dataset %in% c("MSSM", "RUSH", "HBCC", "Aging") ){
			dsArray = dataset
		}else{
			dsArray = c("MSSM", "RUSH", "HBCC")
		}

		for( ds in dsArray){
			write_job( 	df$variable_type[i], 
						ctst_key, 
						ds, 
						df$method[i], 
						df$SampleLevel[i], 
						df$AnnoLevel[i])	
		}
	}
}


