 # Gabriel Hoffman
 # Compute meta-analysis

library(getopt)
# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
  'topTable', 't', 1, "character",
  'code', 'c', 1, "character",
  'method', 'm', 1, "character",
  'outFolder'   , 'o', 1, "character"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

cat("Load libraries...\n")
suppressPackageStartupMessages({
library(openxlsx)
library(tidyverse)
library(arrow)
library(dreamlet)
library(parallel)
})

# library(synapser) synLogin() # synGet("syn53282026")$path 
file = "/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/analysis/contrast_pairs.xlsx"
df_meta = read.xlsx( file )

cat("Read parquet...\n")
df = read_parquet( opt$topTable ) %>%
		select(-AveExpr, -P.Value, -adj.P.Val, -B, -z.std) 

# match code to coefs
#------------------------
coefUniq = unique(df$coef)	

i = which(df_meta$meta == opt$code)

id = c(df_meta$MSSM[i], df_meta$HBCC[i], df_meta$RUSH[i])

# define pattern and coef to subset on 
pattern = paste0("^(", paste0(id[!is.na(id)], collapse="|"), ')' )
coefMatch = grep(pattern, coefUniq, value=TRUE)

grp = c("ID", "assay", "AnnoLevel")

cat("Filtering...\n")
df2 = df %>%
	filter(coef %in% coefMatch) 

cat("Analysis...\n")
res.meta = df2 %>%
		arrange(ID) %>% slice_head(n = 20046) %>%
		meta_analysis(method = opt$method, group=grp) %>%
		mutate( coef = opt$code, Trait = df_meta$Contrast.desc[i])

# write to file
cat("Writing parquet...\n")
outfile = paste0(opt$outFolder, "/topTable_meta_", opt$code, "_", opt$method,".parquet")
write_parquet(res.meta, outfile )

# mashr 
#######

library(mashr)

run_mash_local = function(res, annLvl){

	res = res %>% 
		filter(AnnoLevel == annLvl)

	lvls = unique(res$assay)
	res$assay <- factor(res$assay, lvls)

	# convert to matricies
	B <- dreamlet:::tabToMatrix(res, "estimate")
	S <- dreamlet:::tabToMatrix(res, "std.error")

	# only keep columns with variance in logFC
	cv <- colVars(B, na.rm = TRUE, useNames = TRUE)
	keep <- (cv > 0) & !is.na(cv)
	B <- B[, keep, drop = FALSE]
	S <- S[, keep, drop = FALSE]

	# run mashr on these matricies
	#-----------------------------

	# set up
	# NA's are replaced with beta = 0 with se = 1e6
	data <- mash_set_data(B, S)

	# estimate model parameters
	U.c <- cov_canonical(data)

	# Estimate correlation structure
	V.em <- mash_estimate_corr_em(data, U.c, details = TRUE)

	# copy model to drop NA terms
	model <- V.em$mash.model

	# B has the same ordering as these, so replace corresponding elements with NA
	# this revents non-sensiccal results for coefficients that were originally NA
	idx <- which(is.na(B))
	model$result$PosteriorMean[idx] <- NA
	model$result$PosteriorSD[idx] <- NA
	model$result$NegativeProb[idx] <- NA
	model$result$lfsr[idx] <- NA

	# format results as new object
	res_mash = new("dreamlet_mash_result", list(model = model, logFC.original = B, coefList = opt$code))

	file = paste0(opt$outFolder, "/res_mash_", annLvl, "_", opt$code, '_', opt$method, ".RDS")
	saveRDS(res_mash, file=file)

	ncol(B)
}

lvls = setdiff(unique(res.meta$AnnoLevel), "bulk")

lapply(lvls, function(x){
	message(x)
	run_mash_local(res.meta, x)
	})









