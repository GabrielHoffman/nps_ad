#! /bin/env Rscript

library(getopt)

spec = matrix(c(
  'dataset', 		'd', 1, "character",
  'variable_type', 	'v', 1, "character",
  'AnnoLevel', 		'a', 1, "character",
  'test' , 			't', 1, "character"
), byrow=TRUE, ncol=4)
opt = getopt(spec)


# dataset: MSSM or FULL
# AnnoLevel: class or subclass
params = list(dataset = opt$dataset, #"FULL", 
	variable_type = opt$variable_type, #,"CAT", 
	SampleLevel = "SubID", 
	AnnoLevel = opt$AnnoLevel, #"class", 
	test = opt$test)#, "MSSM_AD__controls")

# Load CONTRASTS and metadata
load("/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/metadata/contrasts_for_dreamlet.Rdata")

suffix2 = with(params, gsub(paste0(dataset, "_", dataset, "_"), paste0(dataset, "_"), paste0(toupper(dataset), "_", test, "_", SampleLevel, "_", AnnoLevel)))

outpath = paste0("/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/analysis/results/grant/", suffix2)

suppressPackageStartupMessages({
library(dreamlet)
library(SingleCellExperiment)
library(tidyverse)
# library(zenith)
library(qvalue)
library(arrow)
})

# read processed data
path = paste0("/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/processAssays/final/")
pattern = paste0(toupper(params$dataset), "_2024-02.*_processAssays_", params$SampleLevel, "_", params$AnnoLevel, ".RDS")
file = dir(path, pattern=pattern, full.names=TRUE)
if( length(file) == 2){                 
  a = file.info(file)$ctime
  file = ifelse(difftime(a[1], a[2]) > 0, file[1], file[2])
}

stopifnot(length(file) ==1)

# read in results from processAssays()
res.proc = readRDS( file )

covariates_base = c("scale(Age)", "Sex", "scale(PMI)", "log(n_genes)", "percent_mito", "mito_genes", "ribo_genes", "mito_ribo", "Brain_bank")
if( params$SampleLevel == "Channel" ){
  covariates_base = c(covariates_base, "(1|SubID)")
}


# combine base covariance with testing variables
covariates = setdiff(union(covariates_base, ctr$covariates_incl), ctr$covariates_excl)

# user-specified formula
if( any(is.null(ctr$contrasts)) || any(is.na(ctr$contrasts)) ){
  form = paste0("~ ", ctr$variable, " + ", paste0(covariates, collapse=" + "))  
  form.vp = form
}else{
  form = paste0("~ 0 + ", ctr$variable, " + ", paste0(covariates, collapse=" + "))
  form.vp = paste0("~ (1|", ctr$variable, ") + ", paste0(covariates, collapse=" + "))
}
form = as.formula(form)
form.vp = as.formula(form.vp)

if( params$SampleLevel == "Channel"){
  # merge subject-level metadata with colData(res.proc)
  metadata_sub = metadata[metadata$SubID %in% colData(res.proc)$SubID,]
  idx = match( colData(res.proc)$SubID, metadata_sub$SubID)
}else{  
  metadata_sub = metadata[metadata$SubID %in% rownames(colData(res.proc)),]
  idx = match( rownames(colData(res.proc)), metadata_sub$SubID)
}

# test ordering
# df = cbind(as.character(colData(res.proc)$SubID), 
#   metadata_sub$SubID[idx])
# all.equal(df[,1], df[,2])

cols = colnames(metadata_sub) %in% colnames(colData(res.proc))
# colData(res.proc) = cbind(colData(res.proc), metadata_sub[idx,!cols])
res.proc@data = cbind(colData(res.proc), metadata_sub[idx,!cols])


if( params$variable_type == "CAT" ){

	# CATEGORICAL
	V = paste(colData(res.proc)[,opt$test], colData(res.proc)$Sex, sep="_")
	V[grep("^NA_", V)] = NA
	tab = table(V)
	cn = paste0("V", names(tab))

	colData(res.proc)$V = V
	form = ~ 0 + V + scale(Age) + scale(PMI) + log(n_genes) + percent_mito +  mito_genes + ribo_genes + mito_ribo 

	# separately test AD in males, females, 
	# then test difference
	contrasts = c(Dx_AD_F= paste(cn[2], "-", cn[4]),
					Dx_AD_M = paste(cn[1], "-", cn[3]),
					Dx_diff = paste("(", paste(cn[2], "-", cn[4]), ") - (", paste(cn[1], "-", cn[3]), ")")	)

	# run dreamlet
	fit = dreamlet( res.proc, form, 
	            contrasts = contrasts)

	coefNames(fit)

	topTable(fit, coef='Dx_diff')

	tabF = topTable(fit, coef='Dx_AD_F', number=Inf) %>% 
			as_tibble %>%
			mutate(se = logFC / t)

	tabM = topTable(fit, coef='Dx_AD_M', number=Inf) %>% 
			as_tibble %>%
			mutate(se = logFC / t)

	df = inner_join(tabF, tabM, by=c("assay", "ID")) %>%
				mutate(z = (logFC.x - logFC.y) / sqrt(se.x^2 + se.y^2),
					p.value = pnorm(z, lower.tail=FALSE),
					FDR = p.adjust(p.value))

	FDR = c()
	for( CT in unique(df$assay) ){
		values = df %>% 
		filter(assay == CT) %>%
		mutate(FDR = p.adjust(p.value)) %>%
		pull(FDR)

		FDR = c(FDR, values)
	}

	df$FDR = FDR
	table(df$FDR < 0.1)

	df$test = opt$test
	df$test = factor(df$test)

	df$AnnoLevel = opt$AnnoLevel
	df$AnnoLevel = factor(df$AnnoLevel)

	file = paste0(outpath, ".parquet")
	write_parquet(df, sink=file)

	# # Plot comparison
	# fig = df %>% 
	# 	ggplot(aes(logFC.x, logFC.y, color=FDR < 0.05)) +
	# 		geom_point() +
	# 		theme_classic() +
	# 		geom_abline(intercept=0, slope=1) + 
	# 		theme(aspect.ratio=1) +
	# 		facet_wrap(~assay) +
	# 		xlab("Disease effect (Female)") +
	# 		ylab("Disease effect (Male)") 
	# ggsave(fig, file="~/www/test.png")


	# # plotVolcano
	# fig = plotVolcano( fit, coef = 'Dx_diff', nGenes=1)
	# ggsave(fig, file="~/www/test.png")

	# go.gs = get_GeneOntology(to="SYMBOL")

	# res = zenith_gsa(fit, go.gs, coefs = "Dx_diff")

	# res2 = res %>%
	# 		as_tibble %>%
	# 		arrange(PValue)

}else{

	##############
	# CONTINUOUS #
	##############

	colData(res.proc)$Trait = colData(res.proc)[,opt$test] %>% as.character %>% as.numeric
	i = !is.na(colData(res.proc)[,opt$test])
	colData(res.proc)$Trait[i] = colData(res.proc)[,opt$test][i] %>% as.character %>% as.numeric

	# MALE
	colData(res.proc)$TraitMale = with(colData(res.proc), Trait * (Sex == "Male"))
	colData(res.proc)$TraitMale[colData(res.proc)$TraitMale==0] = NA

	form = ~ TraitMale + scale(Age) + scale(PMI) + log(n_genes) + percent_mito +  mito_genes + ribo_genes + mito_ribo +Brain_bank
	fitMale = dreamlet( res.proc, form)

	# FEMALE
	colData(res.proc)$TraitFemale = with(colData(res.proc), Trait * (Sex == "Female"))
	colData(res.proc)$TraitFemale[colData(res.proc)$TraitFemale==0] = NA

	form = ~ TraitFemale + scale(Age) + scale(PMI) + log(n_genes) + percent_mito +  mito_genes + ribo_genes + mito_ribo
	fitFemale = dreamlet( res.proc, form)


	tabM = topTable(fitMale, coef='TraitMale', number=Inf) %>% 
			as_tibble %>%
			mutate(se = logFC / t)

	tabF = topTable(fitFemale, coef='TraitFemale', number=Inf) %>% 
			as_tibble %>%
			mutate(se = logFC / t)

	df = inner_join(tabF, tabM, by=c("assay", "ID")) %>%
				mutate(z = (logFC.x - logFC.y) / sqrt(se.x^2 + se.y^2),
					p.value = pnorm(z, lower.tail=FALSE),
					FDR = p.adjust(p.value)) %>%
				arrange(FDR)

	table(df$FDR < 0.1)

	df$test = opt$test
	df$test = factor(df$test)

	df$AnnoLevel = opt$AnnoLevel
	df$AnnoLevel = factor(df$AnnoLevel)

	file = paste0(outpath, ".parquet")
	write_parquet(df, sink=file)

	# df %>%
	# 	head %>%
	# 	data.frame

	# # Plot comparison
	# fig = df %>%
	# 	arrange(-FDR)%>% 
	# 	ggplot(aes(logFC.x, logFC.y, color=FDR < 0.1)) +
	# 		geom_point() +
	# 		theme_classic() +
	# 		geom_abline(intercept=0, slope=1) + 
	# 		theme(aspect.ratio=1) +
	# 		facet_wrap(~assay) +
	# 		xlab("Disease effect (Female)") +
	# 		ylab("Disease effect (Male)") 
	# ggsave(fig, file="~/www/test.png")


	# fig = extractData(res.proc, "EN", 
	# 	cols = c("Sex", "Braak"),
	# 	genes="SYT12") %>%
	# 	filter(!is.na(Braak)) %>%
	# 	ggplot(aes(Braak, SYT12, color=Sex)) +
	# 		geom_point() +
	# 		theme_classic() +
	# 		theme(aspect.ratio=1) +
	# 		geom_smooth(method="lm") +
	# 		scale_color_manual(values=c(Male = "blue", Female="red")) +
	# 		scale_x_continuous(breaks=seq(6))
	    
	# ggsave(fig, file="~/www/test.png")




	# fig = bind_rows(tabM %>% mutate(Sex="Male"), 
	# 	tabF %>% mutate(Sex="Female")) %>%
	# 	inner_join(df %>%
	# 	filter(FDR < 0.05) %>%
	# 	select(assay, ID), 
	# 	by = c("assay", "ID")) %>%
	# 	ggplot(aes(paste(assay, ID), logFC, color=Sex)) +
	# 		geom_hline(yintercept=0, linetype="dashed") +
	# 		geom_point() +
	# 		geom_errorbar(aes(ymin = logFC - 1.96*se, ymax = logFC + 1.98*se), width=0) +
	# 		coord_flip() +
	# 		theme_classic() +
	# 		theme(aspect.ratio=1) +
	# 		scale_color_manual(values=c(Male = "blue", Female="red")) 
	# ggsave(fig, file="~/www/test.png")

}




