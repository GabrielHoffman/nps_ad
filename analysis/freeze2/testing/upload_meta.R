# Gabriel Hoffman
# Jan 16, 2024

# dreamlet results
##################
ml openssl/1.0.2 
R

library(arrow)
library(tidyverse)
library(parallel)
library(dreamlet)

# read data
parent = "/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/analysis/results/"

prefix = "_prs_"
files = dir(parent, pattern="topTable.tsv.gz", recursive=TRUE, full.names=TRUE)
files = grep(prefix, files, value=TRUE)

# files = grep("_prs_", files, value=TRUE, invert=TRUE)

df = mclapply(files, read_tsv, show_col_types=FALSE, mc.cores=12) %>%
		bind_rows

# write results
outfile = paste0(parent, "topTable_combined", prefix, ".tsv.gz")
df %>%
	mutate( logFC = signif(logFC, digits=4),
			AveExpr = signif(AveExpr, digits=4),
			t = signif(t, digits=4),
			P.Value = signif(P.Value, digits=4),
			adj.P.Val = signif(adj.P.Val, digits=4),
			B = signif(B, digits=4),
			z.std = signif(z.std, digits=4)) %>%
	write_tsv(file=outfile )

cmd = paste("ml python; synapse add --parentid syn53144970", outfile)
system(cmd)

outfile2 = paste0(parent, "topTable_combined", prefix, ".parquet")
df %>% 
	write_parquet(outfile2, compression="gzip")	
cmd = paste("ml python; synapse add --parentid syn53144970", outfile2)
system(cmd)


# Joint hypothesis testing
##########################

library(arrow)
library(tidyverse)
library(parallel)
library(dreamlet)

# read data
parent = "/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/analysis/results/"

files = dir(parent, pattern="topTable_joint.tsv.gz", recursive=TRUE, full.names=TRUE)

df = mclapply(files, function(file){
		# joint analysis
		df_joint = read_tsv(file, show_col_types=FALSE)

		# univariate analysis
		file2 = gsub("topTable_joint.tsv.gz", "topTable.tsv.gz", file)
		df = read_tsv(file2, show_col_types=FALSE)

		# annotate df
		df_joint$coef = gsub(".L", "", df$coef[1])
		df_joint$Dataset = df$Dataset[1]
		df_joint$SampleLevel = df$SampleLevel[1]
		df_joint$AnnoLevel = df$AnnoLevel[1]

		df_joint %>%
			select(assay, ID, AveExpr, F, P.Value, adj.P.Val, coef, Dataset, SampleLevel, AnnoLevel)
}, mc.cores=12) %>%
	bind_rows


outfile = paste0(parent, "topTable_combined_joint.tsv.gz")
df %>%
	mutate( AveExpr = signif(AveExpr, digits=4),
			F = signif(F, digits=4),
			P.Value = signif(P.Value, digits=4),
			adj.P.Val = signif(adj.P.Val, digits=4)) %>%
	write_tsv(file=outfile )

cmd = paste("ml python; synapse add --parentid syn53144970", outfile)
system(cmd)

outfile2 = paste0(parent, "topTable_combined_joint.parquet")
df %>% 
	write_parquet(outfile2, compression="gzip")	
cmd = paste("ml python; synapse add --parentid syn53144970", outfile2)
system(cmd)



# mashr analysis
################

library(arrow)
library(tidyverse)
library(parallel)
library(mashr)
library(dreamlet)

# read data
parent = "/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/analysis/results/"

files = dir(parent, pattern="res_mash.RDS", recursive=TRUE, full.names=TRUE)
files = grep("prs_", files, value=TRUE, invert=TRUE)

df_prob = mclapply( files, function(file){

	res = readRDS(file)

	# get probability from lFSR
	prob <- 1 - get_lfsr(res$model)

	info = rev(strsplit(dirname(file), '/')[[1]])

	prob %>%
		as.data.frame %>%
		rownames_to_column(var="ID") %>%
		as_tibble %>%
		pivot_longer(!ID, names_to = "assay", values_to="Probability") %>%	
		mutate(coef = res$coefList, 
			AnnoLevel = info[1], 
			Dataset = info[4]) %>%
		filter(!is.na(Probability))
	}, mc.cores=12)
df_prob = bind_rows(df_prob)

# mashr from meta-analysis
parent = "/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/analysis/results/meta/"
files = dir(parent, pattern="res_mash_.*RDS", recursive=TRUE, full.names=TRUE)
files = grep("prs_", files, value=TRUE, invert=TRUE)

df_prob_meta = mclapply( files, function(file){

	res = readRDS(file)

	# get probability from lFSR
	prob <- 1 - get_lfsr(res$model)

	info = strsplit(basename(file), '_')[[1]]
	info = data.frame(AnnoLevel = info[3], Dataset = "meta")

	prob %>%
		as.data.frame %>%
		rownames_to_column(var="ID") %>%
		as_tibble %>%
		pivot_longer(!ID, names_to = "assay", values_to="Probability") %>%	
		mutate(coef = res$coefList, 
			AnnoLevel = info$AnnoLevel, 
			Dataset = info$Dataset) %>%
		filter(!is.na(Probability))
	}, mc.cores=12)
df_prob_meta = bind_rows(df_prob_meta)

outfile2 = paste0(parent, "mashr_posterior.parquet")
df_prob %>%
	bind_rows(df_prob_meta) %>% 
	write_parquet(outfile2, compression="gzip")	

cmd = paste("ml python; synapse add --parentid syn53144970", outfile2)
system(cmd)



library(arrow)
library(tidyverse)
library(dreamlet)

# df_prob = read_parquet('mashr_posterior.parquet')

# probabilities from first dataset
df1 = df_prob %>%
		filter( Dataset == "MSSM",
			AnnoLevel == "class",
			coef == "c02xAD - c02xControl"	) %>%
		pivot_wider(names_from = "assay", values_from="Probability")

# probabilities from second dataset
df2 = df_prob %>%
		filter( Dataset == "RUSH",
			AnnoLevel == "class",
			coef == "c03xAD - c03xControl"	) %>%
		pivot_wider(names_from = "assay", values_from="Probability")

df_join = inner_join(df1, df2, by="ID") %>%
			select( -coef.x, -AnnoLevel.x, -Dataset.x,
				 	-coef.y, -AnnoLevel.y, -Dataset.y)

# Identify genes differentially expressed in Immune in both datsets, but no other cell types
include = c('Immune.x', 'Immune.y')
exclude = setdiff(colnames(df_join)[-1], include)

res = df_join %>% 
		column_to_rownames(var='ID') %>%
		filter(rowSums(is.na(.)) != ncol(.)) %>%
		compositePosteriorTest(
			include = include,
			exclude = exclude,
			test = "all")

# top genes
head(sort(res, decreasing=TRUE))


# Identify genes differentially expressed in Immune in both datsets, regardless of results in other cell types
include = c('Immune.x', 'Immune.y')
exclude = NULL

res = df_join %>% 
		column_to_rownames(var='ID') %>%
		filter(rowSums(is.na(.)) != ncol(.)) %>%
		compositePosteriorTest(
			include = include,
			exclude = exclude,
			test = "all")

# top genes
head(sort(res, decreasing=TRUE))


# Identify genes differentially expressed in at least one neuron type
include = c('EN.x', 'EN.y', "IN.x", "IN.y")
exclude = setdiff(colnames(df_join)[-1], include)

res = df_join %>% 
		column_to_rownames(var='ID') %>%
		filter(rowSums(is.na(.)) != ncol(.)) %>%
		compositePosteriorTest(
			include = include,
			exclude = exclude)

# top genes
head(sort(res, decreasing=TRUE))





ls /sc/arion/projects/psychAD/NPS-AD/freeze2_rc/analysis/results/HBCC/HBCC_prs_raw_FFM_Agreeableness_GPC.1/SubID/subclass/



# Meta-analysis
###############

cd /sc/arion/projects/psychAD/NPS-AD/freeze2_rc/analysis/

SRC=/sc/arion/projects/CommonMind/hoffman/NPS-AD/work/nps_ad/analysis/freeze2/testing/meta_analysis.R
OUT=/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/analysis/results/meta/
LOG=$(pwd)/results/meta/jobs/../logs/

mkdir -p results/meta/jobs/ $LOG

for N in $(seq -w 01 21)
do
	CODE=$(echo m${N}x)
	for METHOD in $(echo "FE REML RE2C")
	do
echo '#!/bin/bash' > results/meta/jobs/meta_${CODE}_${METHOD}.lsf
echo "#BSUB -J meta_${CODE}_${METHOD}
	#BSUB -P acc_CommonMind
	#BSUB -q premium
	#BSUB -n 1
	#BSUB -R \"span[hosts=1]\"
	#BSUB -R \"rusage[mem=30000]\"
	#BSUB -W 36:00 
	#BSUB -o $LOG/meta_${CODE}_${METHOD}_%J.stdout
	#BSUB -eo $LOG/meta_${CODE}_${METHOD}_%J.stderr
	#BSUB -L /bin/bash
	#BSUB -cwd /sc/arion/projects/psychAD/NPS-AD/freeze2_rc/analysis/

	module purge
	module --ignore-cache load hdf5/1.12.1 libpng/12 R/4.3.3-intel-mkl pandoc/2.6 gsl openssl

	R_LIBS_USER=/hpc/users/hoffmg01/.Rlib/R_433
	R_LIBS=${R_LIBS_USER}:/hpc/packages/minerva-centos7/rpackages/4.3.3-intel-mkl/site-library:/hpc/packages/minerva-centos7/rpackages/bioconductor/3.18

	Rscript $SRC --topTable results/topTable_combined.parquet --code $CODE --method $METHOD --outFolder $OUT" >> results/meta/jobs/meta_${CODE}_${METHOD}.lsf
	done
done


# for PRS
PRS=/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/analysis/prs.labels

for CODE in $(cat $PRS)
do
	for METHOD in $(echo "FE REML RE2C")
	do
echo '#!/bin/bash' > results/meta/jobs/meta_${CODE}_${METHOD}.lsf
echo "#BSUB -J meta_${CODE}_${METHOD}
	#BSUB -P acc_CommonMind
	#BSUB -q premium
	#BSUB -n 1
	#BSUB -R \"span[hosts=1]\"
	#BSUB -R \"rusage[mem=30000]\"
	#BSUB -W 36:00 
	#BSUB -o $LOG/meta_${CODE}_${METHOD}_%J.stdout
	#BSUB -eo $LOG/meta_${CODE}_${METHOD}_%J.stderr
	#BSUB -L /bin/bash
	#BSUB -cwd /sc/arion/projects/psychAD/NPS-AD/freeze2_rc/analysis/

	module purge
	module --ignore-cache load hdf5/1.12.1 libpng/12 R/4.3.3-intel-mkl pandoc/2.6 gsl openssl

	R_LIBS_USER=/hpc/users/hoffmg01/.Rlib/R_433
	R_LIBS=${R_LIBS_USER}:/hpc/packages/minerva-centos7/rpackages/4.3.3-intel-mkl/site-library:/hpc/packages/minerva-centos7/rpackages/bioconductor/3.18

	Rscript $SRC --topTable results/topTable_combined_prs.parquet --code $CODE --method $METHOD --outFolder $OUT" >> results/meta/jobs/meta_${CODE}_${METHOD}.lsf
	done
done

 # results/meta/jobs/*prs*

prs_raw_pd_without_23andMe




ls results/meta/jobs/meta*.lsf | parallel -P1 "bsub < {}"


ls results/meta/jobs/meta*FE*.lsf | parallel -P1 "bsub < {}"

# Combine results
cd /sc/arion/projects/psychAD/NPS-AD/freeze2_rc/analysis/results/meta

library(tidyverse)
library(arrow)

files = dir(".", pattern="topTable_meta_.*_(FE|REML).parquet", full.names=TRUE)

df = lapply(files, read_parquet) %>%
		bind_rows %>%	
		group_by(assay, AnnoLevel) %>%
		mutate(adj.P.Val = p.adjust(p.value))

write_parquet(df, "res_meta.parquet")	
write_tsv(df, "res_meta.tsv.gz")		


system("ml python; synapse add --parentid syn53144970 res_meta.parquet")
system("ml python; synapse add --parentid syn53144970 res_meta.tsv.gz")


# Combine mashr results again
#@@@@@@@@@@@@@@@@@@@@@@@

# crumblr results
#################

library(tidyverse)
library(parallel)

parent = "/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/analysis/results/"

files = dir(parent, pattern="topTable_crumblr.tsv.gz", recursive=TRUE, full.names=TRUE)
files = grep("prs_", files, value=TRUE, invert=TRUE)

df = lapply(files, function(file){

	message(file)

	spl = strsplit(file, '/')[[1]]

	tmp = read_tsv(file, show_col_types=FALSE) 

	if( nrow(tmp) > 0){
		tmp = tmp %>%
		mutate(coef = id, 
			id = NULL,
			Dataset = spl[11],
			SampleLevel = spl[13],
			AnnoLevel = spl[14])
	}else{
		tmp = NULL
	}
	tmp
}) %>%
bind_rows

outfile = paste0(parent, "topTable_combined_crumblr.tsv.gz")
write_tsv(df, file=outfile )

cmd = paste("ml python; synapse add --parentid syn53144970", outfile)

system(cmd)


# Joint hypothesis testing
##########################

library(arrow)
library(tidyverse)
library(parallel)
library(dreamlet)

# read data
parent = "/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/analysis/results/"

files = dir(parent, pattern="topTable_crumblr_joint.tsv.gz", recursive=TRUE, full.names=TRUE)
files = grep("prs_", files, value=TRUE, invert=TRUE)
files = grep("Channel", files, value=TRUE, invert=TRUE)

df = mclapply(files, function(file){
		# joint analysis
		df_joint = read_tsv(file, show_col_types=FALSE)

		# univariate analysis
		file2 = gsub("topTable_crumblr_joint.tsv.gz", "topTable_crumblr.tsv.gz", file)
		df = read_tsv(file2, show_col_types=FALSE)

		s = strsplit(dirname(file),'\\/')[[1]]

		# annotate df
		df_joint$coef = gsub(".L", "", df$id[1])
		df_joint$Dataset = s[11]
		df_joint$SampleLevel = s[13]
		df_joint$AnnoLevel = s[14]

		df_joint %>%
			select(CellType, AveExpr, F, P.Value, adj.P.Val, coef, Dataset, SampleLevel, AnnoLevel)
}, mc.cores=12) %>%
	bind_rows


outfile = paste0(parent, "topTable_crumblr_joint.tsv.gz")
df %>%
	mutate( AveExpr = signif(AveExpr, digits=4),
			F = signif(F, digits=4),
			P.Value = signif(P.Value, digits=4),
			adj.P.Val = signif(adj.P.Val, digits=4)) %>%
	write_tsv(file=outfile )

cmd = paste("ml python; synapse add --parentid syn53144970", outfile)
system(cmd)
