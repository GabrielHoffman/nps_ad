 # Gabriel Hoffman
 # Compute meta-analysis

library(getopt)
# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
  'topTable', 't', 1, "character",
  'code', 'c', 1, "character",
  'outFolder'   , 'o', 1, "character"
), byrow=TRUE, ncol=4)
opt = getopt(spec)


library(openxlsx)
library(tidyverse)
library(arrow)
library(dreamlet)
library(parallel)

# library(synapser) synLogin() # synGet("syn53282026")$path 
file = "/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/analysis/contrast_pairs.xlsx"
df_meta = read.xlsx( file )

df = read_parquet( opt$topTable ) %>%
		select(-AveExpr, -P.Value, -adj.P.Val, -B, -z.std) 

coefUniq = unique(df$coef)	

i = which(df_meta$meta == opt$code)

id = c(df_meta$MSSM[i], df_meta$HBCC[i], df_meta$RUSH[i])

# define pattern and coef to subset on 
pattern = paste0("^(", paste0(id[!is.na(id)], collapse="|"), ')' )
coefMatch = grep(pattern, coefUniq, value=TRUE)

grp = c("ID", "assay", "AnnoLevel")

df2 = df %>%
	filter(coef %in% coefMatch) 

methods = c("FE", "REML", "RE2C")

res = lapply(methods, function(method){

	message(paste(method, opt$code))
	
	# filtering
	df2 %>%
		meta_analysis(method = method, group=grp) %>%
		mutate( coef = x, Trait = df_meta$Contrast.desc[i])
})
names(res) = methods

# write to file
tmp = lapply(names(res), function(method){
	outfile = paste0(opt$outFolder, "/topTable_meta_", opt$code, "_", method,".parquet")
	write_parquet(res[[method]], file=outfile )
	})