

# get latest version
cd /sc/arion/projects/CommonMind/hoffman/NPS-AD/work/nps_ad/analysis/freeze2/
git pull 

1) Merge raw H5AD with metadata to create H5AD files in h5ad_final
`create_freezes/merge_with_metadata.Rmd`

2) Compute pseudobulk, precision weights, etc for each dataset
`preprocess/preprocess.Rmd`

sed 's/$1/RUSH/g' ./submit_preproc.sh | bsub
sed 's/$1/HBCC/g' ./submit_preproc.sh | bsub
sed 's/$1/AGING/g' ./submit_preproc.sh | bsub
sed 's/$1/MSSM/g' ./submit_preproc.sh | bsub
sed 's/$1/FULL/g' ./submit_preproc.sh | bsub

 2a) Write cell type composition to *.tsv files
	`preprocess/write_compositions.R`


# save excludes
head -n 1 exclude_AGING.tsv > exclude.tsv
ls exclude_AGING.tsv exclude_HBCC.tsv exclude_MSSM.tsv exclude_RUSH.tsv | parallel -P1 tail -n +2 >> exclude.tsv

synapse add --parentid syn53144970 exclude.tsv


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
<!-- rm -f jobs/*Channel* -->

# crumblr
######### 
# \ls jobs/RUSH_*crumblr* | parallel -P1 "bsub < {}; sleep .1"
# \ls jobs/HBCC_*crumblr* | parallel -P1 "bsub < {}; sleep .1"
# \ls jobs/MSSM_*crumblr* | parallel -P1 "bsub < {}; sleep .1"

ls jobs/*crumblr* | grep _class | grep SubID | parallel -P1 "bsub < {}; sleep .5"
ls jobs/*crumblr* | grep _subclass | grep SubID | parallel -P1 "bsub < {}; sleep .5"



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


ls jobs/*dreamlet* | grep _bulk | grep SubID | parallel -P1 "bsub < {}; sleep .5"
ls jobs/*dreamlet* | grep _class | grep SubID | parallel -P1 "bsub < {}; sleep .5"
ls jobs/*dreamlet* | grep _subclass | grep SubID | parallel -P1 "bsub < {}; sleep .5"
ls jobs/*dreamlet* | grep _subtype | grep SubID | parallel -P1 "bsub < {}; sleep .5"



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

cat failed_dreamlet.lsf | grep -v -f <( bjobs -w | awk '{print $7}' | grep dreamlet) | parallel -P1 "bsub < {}"


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

7) Combine all results and upload 
cd /sc/arion/projects/psychAD/NPS-AD/freeze2_rc/analysis/

upload_meta.R



covariates_base = c(covariates_base, "percent_mito", "mito_genes", "ribogenes", "mito_ribo")


form = ~ scale(Age) + Sex + scale(PMI) + log(n_genes) + 
    TechPC1 + TechPC2 + TechPC3 + percent_mito + mito_genes + 
    ribo_genes + mito_ribo

dfvp = fitVarPart( res.proc, form, assay=c("EN", "Astro"))

fig = plotVarPart( sortCols(dfvp), label.angle=45)

ggsave(fig, file="~/www/test.png", height=12, width=12)












