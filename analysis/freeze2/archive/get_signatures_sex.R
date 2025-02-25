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
library(tidyverse)
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

	tabF = topTable(fit, coef='Dx_AD_F', number=Inf) %>% 
			as_tibble %>%
			mutate(se = logFC / t)

	tabM = topTable(fit, coef='Dx_AD_M', number=Inf) %>% 
			as_tibble %>%
			mutate(se = logFC / t)

	df = inner_join(tabF, tabM, by=c("assay", "ID")) %>%
				mutate(z = (logFC.x - logFC.y) / sqrt(se.x^2 + se.y^2),
					p.value = pnorm(z, lower.tail=FALSE),
					FDR = p.adjust(p.value)) %>%
				arrange(CT)

	# df_FDR = lapply(unique(df$assay), function(CT){
	# 	df %>% 
	# 	filter(assay == CT) %>%
	# 	mutate(FDR = p.adjust(p.value)) %>%
	# 	select(assay, ID, FDR)
	# })
	# df_FDR = bind_rows(df_FDR)

	# df = inner_join(df, df_FDR, by=c("assay", 'ID'))

	table(df$FDR < 0.1)

	df$test = opt$test
	df$test = factor(df$test)
	df$AnnoLevel = opt$AnnoLevel
	df$AnnoLevel = factor(df$AnnoLevel)

	file = paste0(outpath, ".parquet")
	write_parquet(df, sink=file)

	df2 = topTable(fit, coef='Dx_diff', number=Inf) %>% 
			as_tibble %>%
			mutate(se = logFC / t)

	df2$test = opt$test
	df2$test = factor(df2$test)
	df2$AnnoLevel = opt$AnnoLevel
	df2$AnnoLevel = factor(df2$AnnoLevel)

	file = paste0(outpath, "_diff.parquet")
	write_parquet(df2, sink=file)

	# # Plot comparison
	fig = df %>% 
		ggplot(aes(logFC.x, logFC.y, color=FDR < 0.05)) +
			geom_point() +
			theme_classic() +
			geom_abline(intercept=0, slope=1) + 
			theme(aspect.ratio=1) +
			facet_wrap(~assay) +
			xlab("Disease effect (Female)") +
			ylab("Disease effect (Male)") 
	ggsave(fig, file="~/www/test.png")


	# # plotVolcano
	fig = plotVolcano( fit, coef = 'Dx_diff', nGenes=5, assay="Micro", cutoff = .2)
	ggsave(fig, file="~/www/test.pdf", height=4)

	df %>%
		arrange(p.value) %>%
		head %>%
		select(assay, ID, test, logFC.x, se.x, logFC.y, se.y, z, p.value, FDR) 




	df %>%
		filter(ID == "RAB11B-AS1", assay == "IN_VIP") %>%
		select(assay, ID, test, logFC.x, se.x, logFC.y, se.y, z, p.value, FDR) %>%
		arrange(p.value)


	df %>%
		mutate(FDR = p.adjust(p.value)) %>%
		filter(ID == "LINC01094", assay == "Micro") %>%
		mutate(beta = 2*(logFC.x/se.x^2 - logFC.y/se.y^2)/(1/se.x^2 + 1/se.y^2)) %>%
		mutate(se = sqrt(se.x^2 + se.y^2)) %>%
		data.frame


	df2 %>% 
		filter(ID == "LINC01094", assay == "Micro") 

	fig = df2 %>% 
		group_by(assay) %>%
		summarize(pi1 = 1 - pi0est(P.Value)$pi0) %>%
		mutate(assay = factor(assay, names(cols))) %>%
		ggplot(aes(assay, pi1, fill=assay)) +
			geom_bar(stat="identity") +
			theme_classic() +
			theme(aspect.ratio=1, legend.position="none") +
			scale_y_continuous(limits=c(0, NA), expand=c(0, 0)) + 
			coord_flip() +
			scale_fill_manual(name = "Cell type", values=cols) +
			ylab("Storey's pi1") +
			ggtitle("c02x")
		ggsave(fig, file="~/www/test.pdf", height=10)


	df2 %>% 
		group_by(assay) %>%
		summarize(pi1 = 1 - pi0est(P.Value)$pi0) %>%
		arrange(-pi1)




	fig = extractData(res.proc, "Oligo", 
		cols = c("Sex", "c02x"),
		genes="LILRB3") %>%
		filter(!is.na(c02x)) %>%
		mutate(Disease = as.character(c02x)) %>%
		mutate(Disease = factor(Disease, c("Control", "AD"))) %>%
		ggplot(aes(Disease, `LILRB3`, fill=Disease)) +
			geom_violin() +
			geom_boxplot(width=.1, position = "dodge2", fill="grey80") +
			theme_classic() +
			theme(aspect.ratio=2) +
			geom_smooth(method="lm") +
			scale_fill_manual(values=c(AD = "red3", Control="green3")) +
			facet_wrap(~Sex) +
			ggtitle("LILRB3")
	    
	ggsave(fig, file="~/www/test.pdf")



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




