

# cd /sc/arion/projects/CommonMind/hoffman/sceqtl/sumStats/
# R 

library(tidyverse)
library(UpSetR)
library(arrow)

assay_order = readRDS("assay_order.RDS")

# get colors
#############

# https://rdrr.io/cran/gridpattern/src/R/mean_col.R
#' Compute average color
#'
#' `mean_col()` computes an average color.
#'
#' We currently compute an average color
#' by using the quadratic mean of the colors' RGBA values.
#'
#' @param ... Colors to average
#' @return A color string of 9 characters: `"#"` followed by the
#'         red, blue, green, and alpha values in hexadecimal.
#' @examples
#'  mean_col("black", "white")
#'  mean_col(c("black", "white"))
#'  mean_col("red", "blue")
#' @export
mean_col <- function(...) {
    cols <- unlist(list(...))
	quadratic_mean <- function(x) sqrt(mean(x^2))
    m <- grDevices::col2rgb(cols, alpha=TRUE) / 255.0
    # quadratic mean suggested at https://stackoverflow.com/a/29576746
    v <- apply(m, 1, quadratic_mean)
    grDevices::rgb(v[1], v[2], v[3], v[4])
}

file = "/sc/arion/projects/psychAD/NPS-AD/freeze2_proc/230816_PsychAD_capstone_F1/230921_PsychAD_color_palette.csv"
df_colors = read_csv(file) %>%
				select(category, name, color_hex) %>%
				filter(category %in% c("class", "subclass", "subtype")) %>%
				mutate(Dataset = paste0(category, " / ", name)) %>%
				select(Dataset, color = color_hex) %>%
				bind_rows(tibble(Dataset = "bulk / bulk", color = "grey"))
cols = df_colors$color
names(cols) = df_colors$Dataset
names(cols) = gsub(".* / ", "", names(cols))
cols[["EN_IN"]] = mean_col( cols[["EN"]], cols[["IN"]])


file = c("/sc/arion/projects/CommonMind/zengb02/single_cell_eQTL_for_Gabriel_02212024/eQTL_detection_on_pearsonsubclass_level/cis-eQTL_detection/analysis/evaluate_eQTL_sharing_mashR/get_lfsr/another/lfsr_results_for_psychAD_sceQTL_subclass")

# SUBCLASS
df_subclass = read.table(file, row.names=1, header=TRUE)

df_subclass %>% 
	rownames_to_column("Gene") %>%
	# write_tsv(file="mashr_lfsr_subclass.tsv.gz")


res = apply(df_subclass[,-1], 2, function(x) which(x < 0.05))

ord = assay_order %>%
		grep("^subclass", ., value=TRUE) %>%
		gsub("^subclass / ", '', .)

ord[!ord %in% names(res)]

ord = ord[ord %in% names(res)]

pdf("plots/upset_mashr_subclass.pdf", height=10, width=7)
upset(fromList(res[ord]), order.by = "freq", nsets=6, keep.order=TRUE, sets=ord, nintersects=15)
dev.off()


#
# Strong cell type specificity
##############################
files = c("/sc/arion/projects/CommonMind/zengb02/single_cell_eQTL_for_Gabriel_02212024/eQTL_detection_on_pearsonsubclass_level/cis-eQTL_detection/analysis/evaluate_eQTL_sharing_mashR/get_lfsr/another2/analysis/lfsr_results_for_psychAD_sceQTL_subclass", "/sc/arion/projects/CommonMind/zengb02/single_cell_eQTL_for_Gabriel_02212024/eQTL_detection_on_pearsonsubclass_level/cis-eQTL_detection/analysis/evaluate_eQTL_sharing_mashR/get_lfsr/another2/analysis/lfsr_results_for_psychAD_sceQTL_subclass_additional")

# read in results
df_lfsr = read.table(files[1], row.names=1, header=TRUE) %>%
				rownames_to_column("Gene") %>%
				as_tibble %>%
		mutate(ID = Gene) %>%
		mutate(Gene = gsub("^(\\S+)_.*$", "\\1", Gene))


df2 = read_tsv(files[2]) %>%
	rename(Gene = `cell_colocalized gene` ) %>%
	mutate(ID = gsub("^.* (\\S+)$", "\\1", Gene)) %>%
	mutate(Gene = gsub("^.* (\\S+)_.*$", "\\1", Gene))

df_lfsr = bind_rows(df_lfsr, df2)

# append
df_lfsr %>%
	write_tsv(file="mashr_lfsr_subclass.tsv.gz")



library(dreamlet)

ctlist = list("Astro", "Endo", "Micro", "Oligo", "OPC", "PC", "PVM", "SMC", "VLMC", "EN_L3.*", "EN_L5_6_NP", "EN_L6.*", "EN.*", "IN.*", "IN_LAMP5.*", "IN_PVALB_CHC")

res = lapply(ctlist, function(include){
	message(include)
	df_lfsr2 = df_lfsr %>%
				select(-Gene) %>%
				column_to_rownames("ID") 
	exclude = setdiff(colnames(df_lfsr2), include)

	prob = dreamlet:::.compositePosteriorTest( 1 - df_lfsr2, include = include, exclude = exclude, test="all")

	res = data.frame(prob)
	colnames(res) = paste(include, collapse="_")
	res
})
res = bind_cols(res) %>%
		rownames_to_column("Gene") %>%
		as_tibble

ord.class = assay_order %>%
				grep("^subclass", ., value=TRUE) %>%
				gsub("^subclass / ", "", .) %>%
				c("EN_IN", .)

res_prob = res %>%
			pivot_longer(cols=!Gene, names_to = "subclass", values_to="prob")

res_prob %>%
	arrange(-prob) %>%
	write_tsv(file="specific_eQTL_subclass.tsv.gz")


# Examine two APP tests
# APP_rs128648 is specific to Oligo
# APP_rs2226349 is strongest in Astro, 
#   but has signal in EN_L6_IT_2, IN_LAMP5_RELN, Micro
df_lfsr %>% 
	filter(ID %in% c('APP_rs128648', 'APP_rs2226349')) %>%
	pivot_longer(cols=!Gene&!ID, names_to = "subclass", values_to="lfsr") %>%
	filter(lfsr < 0.01)

res_prob %>% 
	filter(Gene %in% c('APP_rs128648', 'APP_rs2226349')) %>%
	filter(prob > 0.01)



fig = res_prob %>%
	group_by(subclass) %>%
	summarize( nSignif = sum(prob > 0.50)) %>%
	mutate(subclass = factor(subclass, ord.class)) %>%
	ggplot(aes(nSignif, subclass, fill=subclass, label=nSignif)) +
		geom_bar(stat="identity") +
		theme_classic() +
		theme(aspect.ratio=1, legend.position = "none")  +
		scale_x_continuous(limits=c(0, NA), expand=c(0,0)) +
		scale_fill_manual(values = cols) +
		geom_text() +
		xlab("# genes with probability of specific eQTL > 50%")

ggsave(fig, file="plots/specific_eQTL_subclass.pdf")




# Gene heatmap
###############

# Inersect 1) cell type specific genes, and 
# 2) genes that coloc with GWAS


# genes = c("GRID2", "AUTS2", "NXN", "LEPR", "SHANK2", "PDE4B", "LRP8", "NRG3", "PTPRG", "APOE", "TSNARE1", "RERE", "CACNA1A", "CACNA1C", "GRIN2A", "GRIN2B", "BDNF", "SCN2A", "BIN1", "PICALM", "APH1B", "CLU", "STAT3",'INO80D', 
# 'FANCA',
# 'SIRT6',
# 'MAP4K4',
# 'HELLS',
# 'YES1',
# 'IL5',
# 'BRAF')

# genes_shared = c('APH1B', "BRAF", "CLU", "CACNA1C", "IL5", "STAT3")
# genes_specific = c("APOE", "AUTS2", "BIN1", "CACNA1A", "GRIN2B", "EGFR", "PDE4B", "SHANK2", "NXN")
# genes = c(genes_shared, genes_specific, genes)

# genes = res_prob %>%
# 	arrange(-prob) %>%
# 	select(Gene) %>%
# 	distinct %>%
# 	head(400) %>%
# 	pull %>%
# 	c(genes) %>%
# 	unique


# read coloc results
df_coloc = read_tsv("coloc_subclass.tsv", show_col_types=FALSE) %>%
			filter(ppH4 > .5) %>%
			select(Trait, CellType, Gene, ppH4)

# Get set of genes with specific eQTLs at prob > 0 .5
df1 = res_prob %>%
		filter(prob > 0.5) %>%
		dplyr::rename(ID = Gene) %>%
		mutate(Gene = gsub("^(\\S+)_.*", "\\1", ID)) %>%
		select(Gene, ID, subclass) %>%
		dplyr::rename(CellType.specific = subclass)


# Get gene cell type pairs for 
# 1) Genes that coloc in any cell type and 
# 2) have a specific eQTL in any cell type
df_merge = df_coloc %>%
	# select(Gene, CellType) %>%
	dplyr::rename(CellType.coloc = CellType) %>%
	distinct %>%
	inner_join(df1, relationship = "many-to-many") %>%
	filter(CellType.coloc != CellType.specific) %>%
	filter(!grepl("\\.\\*", CellType.specific)) %>%
	arrange(Gene) 

df_merge %>%
	write_tsv(file="merge_coloc_mashr.tsv")


df1 %>%
	filter(Gene == 'PFKFB2')

df_coloc %>%
	filter(Gene == 'PFKFB2')


df1 %>%
	filter(Gene == 'APP')

df_coloc %>%
	filter(Gene == 'APP')

df_merge %>%
	filter(Gene == 'APP')


# cell type specific eQTLs that also coloc with GWAS
res_spec_coloc = res_prob %>%
		dplyr::rename(ID = Gene) %>%
		mutate(Gene = gsub("^(\\S+)_.*", "\\1", ID)) %>%
		mutate(CellType = factor(subclass, c(ord.class))) %>%
		select(-subclass) %>%
		full_join(df_coloc, by=c("Gene"), relationship = "many-to-many") %>%
		filter(!is.na(Trait)) %>%
		# filter(prob > 0.5) %>%
		arrange(Gene) 

res_spec_coloc %>% 
	filter(CellType.x == CellType.y) %>%
	filter(prob > 0.5) %>%
	write_tsv(file="coloc_and_subclass_specific.tsv")



res_spec_coloc  %>% 
	filter(Gene == "APP")

# Plot lfsr examples
#####################

library(aplot)

df_lfsr = read_tsv("mashr_lfsr_subclass.tsv.gz")


df_shared = res_spec_coloc %>%
			filter(prob > 0.5) %>% 
			filter(CellType.x == CellType.y)


res_spec_coloc %>% 
	filter(Gene == "RERE")


df_coloc %>%
	filter(Gene == "AUTS2") %>%
	arrange(-ppH4)




# select genes
genes1 = df_shared %>%
			pull(Gene)

# genes = res_prob %>%
# 		dplyr::rename(ID = Gene) %>%
# 		mutate(Gene = gsub("^(\\S+)_.*", "\\1", ID)) %>%
# 		mutate(CellType = factor(subclass, c(ord.class)))  %>%
# 		arrange(-prob) %>%
# 		head(20) %>%
# 		pull(Gene)

ord.class = assay_order %>%
				grep("^subclass", ., value=TRUE) %>%
				gsub("^subclass / ", "", .) 
		


genes = c(genes1, 'CACNA1C', 'DOCK1', 'GRIA1', 'INPP5D', 'NALCN', 'SERPINB1', 'PSD3', 'WNT5B', "EGFR", "TLE4", "BDNF", "DRD2", "GRIN2A", "RERE", "PTPRG",  "CLU", "SHANK2", "APP")


M = df_lfsr %>%
	filter(Gene %in% genes) %>%
	mutate(Gene = ID) %>%
	select(-ID) %>%
	column_to_rownames('Gene')

hcl = hclust(dist(-log10(M+1e-4)))	
id.ord = hcl$labels[hcl$order]

id.ord = c(
"BIN1_rs6733839",
"ATXN1_rs56374431",
"MAX_rs1957949",
"EPHA1-AS1_rs11765305",
"INPP5D_rs10933431",
"INPP5D_rs4663834",
"SERPINB1_rs1350850",
"DOCK1_rs928822631",
"GALNT6_rs3782473",
"SH3TC2_rs28173",
"TLE4_rs17082511",
"TLE4_rs2807303",
"PSD3_rs2638657",
"PSD3_rs2012342",
"SLC4A8_rs3782473",
"SHANK2_rs1660884",
"APP_rs128648",
"APP_rs2226349",
"MYT1L_rs3749051",
"NALCN_rs61973700",
"RERE_rs301796",
"GRIA1_rs140764433",
"NALCN_rs4329782",
"NALCN_rs113644271",
"EGFR_rs74504435",
"CUL3_rs10933068",
"WNT5B_rs7969361",
"WNT5B_rs2270037",
"PFKFB2_rs11585761",
"GRM1_rs2300627",
"AIG1_rs142387624",
"CXCL14_rs2214272",
"CAMK1D_rs10906209",
"PTPRG_rs638552",
"BDNF_rs12417583",
"DRD2_rs10736470",
"CACNA1C_rs2429147",
"DRD2_rs2089652",
"GRIN2A_rs12932206",
"CLU_rs894020")

# hcl$labels[hcl$order][!hcl$labels[hcl$order] %in% id.ord]

exclude = c("NALCN_rs61973700","TLE4_rs4877145", "MYT1L_rs6727254")

df_plot = df_lfsr %>%
	filter(Gene %in% genes) %>%
	pivot_longer(cols=!Gene&!ID, names_to = "class", values_to="lfsr") %>%
	mutate(SNP = gsub("^(\\S+)_(\\S+)$", "\\2", ID)) %>%
	mutate(class = factor(class, ord.class)) %>%
	filter( ! ID %in% exclude, lfsr < .98) %>%
	mutate(ID = factor(ID, rev(id.ord)))

df_plot %>%
	filter(is.na(ID))

fig1 = df_plot %>% 
	ggplot(aes(class, ID,  color = -log10(lfsr+1e-4), size = -log10(lfsr+1e-4), label=ifelse(lfsr < 0.05, 'x', ''))) +
		geom_point() +
		theme_classic() + 
		coord_equal() +
		scale_size_area() +
		scale_color_gradient(low="white", high="red", limits=c(0, 4)) +
		scale_x_discrete(guide = guide_axis(angle = 90)) +
		geom_text(color="black", vjust=0.5, hjust=0.5) 


fig2 = tibble(ID = id.ord) %>%
	left_join(res_spec_coloc)  %>%
	mutate(ID = factor(ID, rev(id.ord)), y = 1) %>%
	filter( ! ID %in% exclude) %>%
	select(Trait, ID) %>%
	distinct %>%
	ggplot(aes( Trait, ID, fill=Trait)) +
		geom_tile() +
		coord_equal() +
		scale_fill_brewer(palette="Set1")

fig = fig1 + fig2


ggsave(fig, file="plots/gene_heatmap.pdf", height=10, width=14)














# 

genes = c('NALCN', 'WNT5B', 'DOCK1', "GRIA1", 'APP', 'CACNA1C', 'PSD3', "INPP5D", "SERPINB1", "TLE4", res_spec_coloc$Gene)
genes = sort(unique(genes))

df_sub = df_join %>%	
	filter(Gene %in% genes) %>%
	mutate(Gene = factor(Gene, genes)) %>%
	mutate(CellType = factor(CellType, c(ord.other, ord.subclass))) %>%
	mutate(variant = gsub("^\\S*_(\\S+)$", "\\1", ID))


files3 = dir('/sc/arion/projects/CommonMind/hoffman/sceqtl/sumStats/', pattern="meta-.*-subclass.*parquet", full.names=TRUE)

dSet = open_dataset(files3)

# get target SNP
df_merge = dSet %>%
		rename(Gene = gene) %>%
		filter(Gene %in% genes)  %>%
		inner_join(df_sub, by=c("Gene", "CellType", "variant")) %>%
			collect

# get target beta
df_target = df_merge %>%
	# filter(Gene == "APP") %>% 
	group_by(ID) %>%
	summarize(beta.target = beta[which.max(abs(z))])


fig = df_merge %>%
	select(-Source, -Level, -chr, -start, -end, -nonassessed_allele, -assessed_allele, -isTopHit, -se, -Trait, -ppH4 ) %>%
	inner_join(df_target)  %>%
	filter(Gene == "APP") %>% 
	# filter(variant == "rs2226349") %>%
	# filter(variant == "rs128648") %>%
	# filter(CellType %in% c("Oligo", "Astro", "IN_LAMP5_RELN")) %>%
	mutate(CellType = factor(CellType, c(ord.other, ord.subclass))) %>%
	mutate(score =if_else(sign(beta) == sign(beta.target), -log10(lfsr + 1e-4), -log10(pmax(1-lfsr- 1e-4, lfsr+ 1e-4)) )) %>%
	ggplot(aes(ID, CellType, color = score, size=score, label=ifelse(lfsr < 0.05, "x", ''))) +
		geom_point() +
		theme_classic() + 
		coord_equal() +
		scale_size_area() +
		scale_color_gradient(low="white", high="red", limits=c(0, 4)) +
		scale_x_discrete(guide = guide_axis(angle = 90)) +
		geom_text(color="black", hjust=0.5, vjust=0.5) 

ggsave(fig, file="plots/gene_heatmap.pdf", width=20)




df_sub %>%
	filter(!is.na(ppH4))


# DRD2 mashr in Astro, IN_VIP, EN_L3_5_IT_3
# coloc with SCZ in EN_L2_3_IT, EN_L3_5_IT_2




df_subclass %>% 
	filter(grepl("^APP_", ID)) %>% 
	data.frame





cd /sc/arion/projects/CommonMind/hoffman/sceqtl/sumStats/
R 

suppressPackageStartupMessages({
library(tidyverse)
library(arrow)
library(ggplot2)
})


plotGenePanels = function(files, GENE, CT, SNP, window = 2e6, min.p=1e-7){

	dSet = open_dataset(files)

	# Forrest plot
	#--------------

	# get target SNP
	df_sub = dSet %>%
				filter(
					# isTopHit, 
					variant == SNP,
					CellType == CT, 
					gene == GENE) %>%
				collect

	df = dSet %>%
			filter(gene == GENE, 
				variant == df_sub$variant) %>%
				collect


	ord.class = assay_order %>%	
				grep("^subclass", ., value=TRUE) %>%
				gsub("^subclass / ", "", .) 

	fig.forrest = df %>%
		mutate(CellType = factor(CellType, ord.class)) %>%
		ggplot(aes(beta, CellType, color = CellType)) +
			geom_point() +
			geom_errorbar(aes(xmin = beta - 1.96*se, xmax = beta + 1.96*se), width=0) +
			theme_classic() +
			theme(legend.position="none", plot.title = element_text(hjust = 0.5)) +
			scale_color_manual(values = cols) +
			geom_vline(xintercept=0) +
			ggtitle(paste(GENE, df_sub$variant, sep=' - '))

	# zoom manhattan plot
	#--------------------

	df_sub = dSet %>%
				filter(gene == GENE) %>%
				collect %>%
				mutate(p.value = 2*pnorm(abs(z), lower.tail=FALSE))

	# filter by window width
	center = mean(df_sub$start)
	df_sub = df_sub %>%
		filter(start > center - window, start < center + window) %>%
		arrange(-p.value)

	CTs = df_sub %>%
		group_by(CellType) %>%
		summarize(p.min = min(p.value)) %>%
		filter(p.min < min.p) %>%
		pull(CellType)

	fig.mht = df_sub %>%
			filter( CellType %in% CTs) %>%
			mutate(isTarget = variant == SNP) %>%
			arrange(isTarget) %>%
			ggplot(aes(start, -log10(p.value), color=isTarget)) +
				geom_point() +
				theme_classic() +
				theme(plot.title = element_text(hjust = 0.5), legend.position='none') +
				facet_wrap(~CellType, ncol=1, scales="free_y") +
				ggtitle(paste(GENE, '/', SNP)) +
				scale_color_manual(values=c("black", "red"))

	# Expression level
	#-----------------

	fig.expr = read_tsv("expressionSpecificity_subclass.tsv.gz", show_col_types=FALSE, progress=FALSE) %>%
		filter(Gene == GENE) %>%
		pivot_longer(!Gene &!totalCPM) %>%
		mutate(CPM = totalCPM * value) %>%		
		mutate(name = factor(name, ord.class)) %>%
		filter(name %in% df$CellType) %>%
		ggplot(aes(Gene, name, fill=log2(CPM))) +
			geom_tile() +
			theme_classic() + 
			coord_equal() +
			theme(plot.title = element_text(hjust = 0.5)) +
			scale_fill_gradient(low="white", high="red", limits=c(0, NA)) +
			xlab('') +
			ylab('') +
			ggtitle(GENE)

	list(fig.forrest = fig.forrest, 
		fig.mht = fig.mht,
		fig.expr = fig.expr)
}			


df_lfsr = read_tsv("mashr_lfsr_subclass.tsv.gz")
df_overlap = read_tsv("specific_and_coloc.tsv")


# Get set of top variants
#-------------------------
# files1 = dir('/sc/arion/projects/CommonMind/hoffman/sceqtl/sumStats/', pattern="-bulk", full.names=TRUE)
# files2 = dir('/sc/arion/projects/CommonMind/hoffman/sceqtl/sumStats/', pattern="-class", full.names=TRUE)
files3 = dir('/sc/arion/projects/CommonMind/hoffman/sceqtl/sumStats/', pattern="meta-.*-subclass.*parquet", full.names=TRUE)

lst = list(,
		c(GENE = "CACNA1A", CT = "Micro"),
		c(GENE = "APOE", CT = "Micro"),
		c(GENE = "BDNF", CT = "EN_L2_3_IT"),
		c(GENE = "CLU", CT = "Micro"),
		c(GENE = "PDE4B", CT = "Oligo"),
		c(GENE = "SHANK2", CT = "Oligo"),
		c(GENE = "NXN", CT = "IN_PVALB_CHC"),
		c(GENE = "STAT3", CT = "IN_PVALB_CHC"),
		c(GENE = "EGFR", CT = "Astro", SNP="rs74504435"),
		c(GENE = "RERE", CT = "Astro"),
		c(GENE = "APP", CT = "Oligo", SNP='rs128648'),	
		c(GENE = "APP", CT = "Astro", SNP='rs2226349'),		
		c(GENE = "DRD2", CT = "IN_VIP"),
		x)

lst2 = lapply(seq(nrow(df_overlap)), function(i){
	c(GENE = df_overlap$Gene[i], GENE = df_overlap$CellType[i])
	})

for(x in  append(lst, lst2)){
	GENE = x[1]
	CT = x[2]
	SNP = x[3]
	message(GENE)

	res = plotGenePanels(files3, GENE, CT, SNP, 2e6, min.p=1e-5)

	file = paste0("plots/details/", GENE, "_", CT, "_", SNP, "_forrest.pdf")
	ggsave(res$fig.forrest, file=file)
	file = paste0("plots/details/", GENE, "_", CT, "_", SNP, "_manhattan.pdf")
	ggsave(res$fig.mht, file=file, height=10, width=4)

	file = paste0("plots/details/", GENE, "_", CT, "_", SNP, "_expression.pdf")
	fig = aplot::plot_list(res$fig.forrest + theme(aspect.ratio=1), res$fig.expr )
	ggsave(fig, file=file, height=7, width=9)
}


# Mashr for Bryois and De Jager
###############################

library(tidyverse)

# file = "/Users/gabrielhoffman/Library/CloudStorage/Dropbox/projects/PsychAD3.0/supplementary_figures/mashR_for_Bryois_and_DeJager/Bryois/lfsr_results_for_Bryois_eQTL.tsv"
# df_bryois = read.table( file, sep="\t", row.names=1 )

# file = "/Users/gabrielhoffman/Library/CloudStorage/Dropbox/projects/PsychAD3.0/supplementary_figures/mashR_for_Bryois_and_DeJager/DeJager/lfsr_results_for_DeJager_eQTL.tsv"


files = c("/sc/arion/projects/CommonMind/zengb02/eQTL_catalogue/down_eQTL_summary_results/brain_sceQTL_DeJager/evaluate_eQTL_sharing_mashR/get_lfsr/another/analysis/lfsr_results_for_DeJager_sceQTL_subclass", "/sc/arion/projects/CommonMind/zengb02/eQTL_catalogue/down_eQTL_summary_results/brain_sceQTL_DeJager/evaluate_eQTL_sharing_mashR/get_lfsr/another/analysis/lfsr_results_for_DeJager_sceQTL_subclass_additional")
df_dj = lapply(files, function(file){
	read.table( file, sep="\t", row.names=1, header=TRUE )
})
rownames(df_dj[[2]]) = gsub("^(\\S+) ", "", rownames(df_dj[[2]]))
df_dj = bind_rows(df_dj) %>%
			rownames_to_column("ID") %>%
			as_tibble %>%
			mutate(Gene = gsub("^(\\S+)_(\\S+)_(\\S+)$", "\\1", ID))%>%
			mutate(SNP = gsub("^(\\S+)_(\\S+)_(\\S+)$", "\\3", ID))



df_dj %>%
	filter(Gene == "APP") %>%
	data.frame




df_dj %>%
	filter(Gene == "ADARB1")


df_dj %>%
	filter(Gene == "EGFR")


df_dj %>%
	filter(Gene == "GRM1")










df_merge = inner_join(df_class, df_dj, by="Gene")

df_merge %>%
	with(cor(-log10(EN+1e-16), -log10(Exc+1e-16), method="sp"))




fig = df_merge %>%
	ggplot(aes(-log10(EN+1e-6), -log10(Exc+1e-6))) +
		geom_point() +
		theme_classic() +
		theme(aspect.ratio=1) +
		coord_equal(xlim=c(0, 6.2)), ylim=c(0, 6.2), expand=FALSE)

ggsave(fig, file="mashr_concordance.pdf")



fig = df_merge %>%
	ggplot(aes(1 - Astro, 1 - Exc)) +
		geom_point() +
		theme_classic() +
		theme(aspect.ratio=1) +
		coord_equal(xlim=c(0, 1), ylim=c(0, 1), expand=FALSE)

ggsave(fig, file="mashr_concordance.pdf")


cutoffs = seq(1e-4, 1-1e-4, length.out=1000)

res = lapply(cutoffs, function(cutoff){
	i = which(df_merge$Astr < cutoff)
	j = df_merge$Exc[i] < cutoff

	tab = with(df_merge, table(EN < cutoff, Exc < cutoff))

	count = sum(j)
	total = length(i)
	tibble(cutoff, count, total, OR = fisher.test(tab)$estimate)
	})
res = bind_rows(res)


fig = res %>%
		ggplot(aes(cutoff, count / total)) +
		geom_point() +
		theme_classic() +
		theme(aspect.ratio=1) +
		scale_y_continuous(limits=c(0,1))


ggsave(fig, file="mashr_concordance.pdf")



fig = res %>%
		ggplot(aes(cutoff, log2(OR))) +
		geom_point() +
		theme_classic() +
		theme(aspect.ratio=1) 

ggsave(fig, file="mashr_concordance.pdf")













