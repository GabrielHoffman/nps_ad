

# get latest version
cd /sc/arion/projects/CommonMind/hoffman/NPS-AD/work/nps_ad/analysis/freeze2/
ml python git pandoc 
ml gcc/11.2.0
git pull 

1) Merge raw H5AD with metadata to create H5AD files in h5ad_final
`create_freezes/merge_with_metadata.Rmd`

2) Compute pseudobulk, precision weights, etc for each dataset
`preprocess/preprocess.Rmd`

 2a) Write cell type composition to *.tsv files
	`preprocess/write_compositions.R`

sed 's/$1/RUSH/g' ./submit_preproc.sh | bsub
sed 's/$1/HBCC/g' ./submit_preproc.sh | bsub
sed 's/$1/MSSM/g' ./submit_preproc.sh | bsub
sed 's/$1/AGING/g' ./submit_preproc.sh | bsub

3) Hypothesis testing for each variable
`analysis/freeze2/testing/write_contrast_jobs.R`
run code is in `run_dreamlet_contrasts.Rmd`

4) Create jobs 
# path for analysis code
cd /sc/arion/projects/CommonMind/hoffman/NPS-AD/work/nps_ad

# Create jobs
R
source("/sc/arion/projects/CommonMind/hoffman/NPS-AD/work/nps_ad/analysis/freeze2/testing/write_contrast_jobs.R")

5) Run jobs
cd /sc/arion/projects/psychAD/NPS-AD/freeze2_rc/analysis/

# crumblr
######### 
# \ls jobs/RUSH_*crumblr* | parallel -P1 "bsub < {}; sleep 1"
# \ls jobs/HBCC_*crumblr* | parallel -P1 "bsub < {}; sleep 1"
# \ls jobs/MSSM_*crumblr* | parallel -P1 "bsub < {}; sleep 1"

# Need more work on ordinal variables:
RUSH_BRAAK_SubID_subtype_crumblr

should this be continuous as well?
treeTest for multiple coef

# dreamlet
##########

# bsub < jobs/MSSM_AD__controls_SubID_subclass_dreamlet.lsf

# \ls jobs/RUSH_*dreamlet* | grep SubID | parallel -P1 "bsub < {}; sleep 1"
# \ls jobs/HBCC_*dreamlet* | grep SubID | parallel -P1 "bsub < {}; sleep 1"
# \ls jobs/MSSM_*dreamlet* | grep SubID | parallel -P1 "bsub < {}; sleep 1"

# Note that these are duplicates, only one is evaluated
# the other ends up empty and fails
# results/RUSH/RUSH_AD__controls/SubID/class
# results/RUSH/AD__controls/SubID/class


isEmptyFolder () {
	N=$(ls -A $1 | wc -l)
	if [ $ > 0 ]
}
export -f isEmptyFolder

# results in 
# cd /sc/arion/projects/psychAD/NPS-AD/freeze2_rc/analysis/

# list of dreamlet jobs that did not produce topTable.tsv.gz
comm -13 <(ls results/*/*/SubID/*/topTable.tsv.gz | parallel -P1 dirname | sort -u) <(ls results/*/*/SubID/* | grep ':' | sed 's/://g' | sort -u) | sed 's/results\///g' | tr '/' '_' | sed 's/MSSM_MSSM/MSSM/g' | sed 's/RUSH_RUSH/RUSH/g'  | sed 's/HBCC_HBCC/HBCC/g' | parallel -P1 echo "jobs/{}_dreamlet.lsf" > failed_dreamlet.lsf

wc failed_dreamlet.lsf

cat failed_dreamlet.lsf | parallel -P1 "bsub < {}"


# list of crumblr jobs that did not produce topTable.tsv.gz
comm -13 <(ls results/*/*/SubID/*/topTable_crumblr.tsv.gz | parallel -P1 dirname | sort -u) <(ls results/*/*/SubID/* | grep ':' | sed 's/://g' | sort -u) | sed 's/results\///g' | tr '/' '_' | sed 's/MSSM_MSSM/MSSM/g' | sed 's/RUSH_RUSH/RUSH/g' | sed 's/HBCC_HBCC/HBCC/g' | parallel -P1 echo "jobs/{}_crumblr.lsf" | grep -v bulk > failed_crumblr.lsf


wc failed_crumblr.lsf

cat failed_crumblr.lsf | parallel -P1 "bsub < {}"

ls results/*/*/SubID/*/topTable_crumblr.tsv.gz | parallel -P1 dirname | sort -u | grep HBCC | grep SCZ_CTRL
ls results/*/*/SubID/* | grep ':' | sed 's/://g' | sort -u | grep HBCC | grep SCZ_CTRL



SCZ_CTRL_SubID_class_crumblr


ls results/*/*/SubID/* | grep ':' | sed 's/://g' | sort -u | parallel -P1 isEmptyFolder {}

6) HTML results
https://hoffmg01.hpc.mssm.edu/PsychAD_analysis/


c("scale(Age)", "Sex", "scale(PMI)", "log(n_genes)", "(1|BatchID)") 

scale(Age) + Sex + scale(PMI) + log(n_genes) + (1|BatchID)

7) Combine all results

# dreamlet results
##################

library(tidyverse)
library(parallel)

parent = "/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/analysis/results/"

files = dir(parent, pattern="topTable.tsv.gz", recursive=TRUE, full.names=TRUE)

df = mclapply(files, read_tsv, show_col_types=FALSE, mc.cores=12) %>%
		bind_rows

outfile = "topTable_combined.tsv.gz"

df %>%
	mutate( logFC = signif(logFC, digits=4),
			AveExpr = signif(AveExpr, digits=4),
			t = signif(t, digits=4),
			P.Value = signif(P.Value, digits=4),
			adj.P.Val = signif(adj.P.Val, digits=4),
			B = signif(B, digits=4),
			z.std = signif(z.std, digits=4)) %>%
	write_tsv(file=outfile )

system("ml python; synapse add --parentid syn53144970 topTable_combined.tsv.gz")

# crumblr results
#################

library(tidyverse)
library(parallel)

parent = "/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/analysis/results/"

files = dir(parent, pattern="topTable_crumblr.tsv.gz", recursive=TRUE, full.names=TRUE)

df = lapply(files, function(file){

	spl = strsplit(file, '/')[[1]]

	read_tsv(file, show_col_types=FALSE) %>%
		mutate(coef = id, 
			id = NULL,
			Dataset = spl[11],
			SampleLevel = spl[13],
			AnnoLevel = spl[14])
}) %>%
bind_rows

outfile = "topTable_combined_crumblr.tsv.gz"

write_tsv(df, file=outfile )

system("ml python; synapse add --parentid syn53144970 topTable_combined_crumblr.tsv.gz")








