

cd /sc/arion/projects/CommonMind/hoffman/sceqtl/sumStats/
R 

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


# files = system("ls /sc/arion/projects/CommonMind/zengb02/single_cell_eQTL_for_Gabriel_02212024/eQTL_detection_on_pearsonsubclass_level/cis-eQTL_detection/analysis/evaluate_eQTL_sharing_mashR/get_lfsr/lfsr_results_for_psychAD_sceQTL_*", intern=TRUE)

files = c("/sc/arion/projects/CommonMind/zengb02/single_cell_eQTL_for_Gabriel_02212024/eQTL_detection_on_pearsonclass_level/cis-eQTL_detection/analysis/evaluate_eQTL_sharing_mashR/get_lfsr/another2/lfsr_results_for_psychAD_sceQTL_class", "/sc/arion/projects/CommonMind/zengb02/single_cell_eQTL_for_Gabriel_02212024/eQTL_detection_on_pearsonsubclass_level/cis-eQTL_detection/analysis/evaluate_eQTL_sharing_mashR/get_lfsr/another/lfsr_results_for_psychAD_sceQTL_subclass")

# CLASS
df_class = read.table(files[1], row.names=1, header=TRUE) %>%
				rownames_to_column("Gene") %>%
				as_tibble
colnames(df_class) = colnames(df_class) %>% 
						gsub("_neurons", '', .) %>% 
						gsub("EX", "EN", .) %>% 
						gsub("Microglia", "Immune", .) %>% 
						gsub("Astrocytes", "Astro", .) %>% 
						gsub("Oligodendrocyte", "Oligo", .) 


res = apply(df_class[,-1], 2, function(x) which(x < 0.05))

ord = assay_order %>%
		grep("^class", ., value=TRUE) %>%
		gsub("^class / ", '', .)

# ord = c("IN", "EN", "Oligo", "OPC", "Astro", "Immune")
ord = assay_order %>%
		grep("^class", .,value=TRUE) %>%
		gsub("class / ", "", .)

pdf("plots/upset_class.pdf", height=5, width=7)
upset(fromList(res[ord]), order.by = "freq", nsets=8, keep.order=TRUE, sets=ord, nintersects=15)
dev.off()


# SUBCLASS
# files = c("/sc/arion/projects/CommonMind/zengb02/single_cell_eQTL_for_Gabriel_02212024/eQTL_detection_on_pearsonsubclass_level/cis-eQTL_detection/analysis/evaluate_eQTL_sharing_mashR/get_lfsr/another2/analysis/lfsr_results_for_psychAD_sceQTL_subclass", "/sc/arion/projects/CommonMind/zengb02/single_cell_eQTL_for_Gabriel_02212024/eQTL_detection_on_pearsonsubclass_level/cis-eQTL_detection/analysis/evaluate_eQTL_sharing_mashR/get_lfsr/another2/analysis/lfsr_results_for_psychAD_sceQTL_subclass_additional")

df_subclass = read.table(files[2], row.names=1, header=TRUE)

res = apply(df_subclass[,-1], 2, function(x) which(x < 0.05))

ord = assay_order %>%
		grep("^subclass", ., value=TRUE) %>%
		gsub("^subclass / ", '', .)

ord[!ord %in% names(res)]

ord = ord[ord %in% names(res)]

pdf("plots/upset_subclass.pdf", height=10, width=7)
upset(fromList(res[ord]), order.by = "freq", nsets=6, keep.order=TRUE, sets=ord, nintersects=15)
dev.off()


#
# Strong cell type specificity
##############################
files = c("/sc/arion/projects/CommonMind/zengb02/single_cell_eQTL_for_Gabriel_02212024/eQTL_detection_on_pearsonsubclass_level/cis-eQTL_detection/analysis/evaluate_eQTL_sharing_mashR/get_lfsr/another2/analysis/lfsr_results_for_psychAD_sceQTL_subclass", "/sc/arion/projects/CommonMind/zengb02/single_cell_eQTL_for_Gabriel_02212024/eQTL_detection_on_pearsonsubclass_level/cis-eQTL_detection/analysis/evaluate_eQTL_sharing_mashR/get_lfsr/another2/analysis/lfsr_results_for_psychAD_sceQTL_subclass_additional")

df_subclass = read.table(files[1], row.names=1, header=TRUE) %>%
				rownames_to_column("Gene") %>%
				as_tibble %>%
		mutate(ID = Gene) %>%
		mutate(Gene = gsub("^(\\S+)_.*$", "\\1", Gene))

df2 = read_tsv(files[2]) %>%
	rename(Gene = `cell_colocalized gene` ) %>%
	mutate(ID = gsub("^.* (\\S+)$", "\\1", Gene)) %>%
	mutate(Gene = gsub("^.* (\\S+)_.*$", "\\1", Gene))

df_subclass = bind_rows(df_subclass, df2)

df_lfsr = df_class %>%
			column_to_rownames("Gene")

ctlist = as.list(colnames(df_lfsr))
ctlist[[length(ctlist)]] = c("EN", "IN")

res = lapply(ctlist, function(include){

	exclude = setdiff(colnames(df_lfsr), include)

	prob = dreamlet:::.compositePosteriorTest( 1 - df_lfsr, include = include, exclude = exclude, test="all")

	res = data.frame(prob)
	colnames(res) = paste(include, collapse="_")
	res
})
res = bind_cols(res) %>%
		rownames_to_column("Gene") %>%
		as_tibble

ord.class = assay_order %>%
				grep("^class", ., value=TRUE) %>%
				gsub("^class / ", "", .) %>%
				c("EN_IN", .)

res_prob = res %>%
			pivot_longer(cols=!Gene, names_to = "class", values_to="prob")

res_prob %>%
	arrange(-prob) %>%
	write_tsv(file="specific_eQTL_class.tsv.gz")

fig = res_prob %>%
	group_by(class) %>%
	summarize( nSignif = sum(prob > 0.50)) %>%
	mutate(class = factor(class, ord.class)) %>%
	ggplot(aes(nSignif, class, fill=class, label=nSignif)) +
		geom_bar(stat="identity") +
		theme_classic() +
		theme(aspect.ratio=1, legend.position = "none")  +
		scale_x_continuous(limits=c(0, NA), expand=c(0,0)) +
		scale_fill_manual(values = cols) +
		geom_text() +
		xlab("# genes with probability of specific eQTL > 50%")

ggsave(fig, file="plots/specific_eQTL_class.pdf")


# SUBCLASS
#-----------------------

df_lfsr = df_subclass %>%
			select(!Gene)  %>%
			column_to_rownames('ID') 

ctlist = list("Astro", "Endo", "Micro", "Oligo", "OPC", "PC", "PVM", "SMC", "VLMC", "EN_L3.*", "EN_L5_6_NP", "EN_L6.*", "EN.*", "IN.*", "IN_LAMP5.*", "IN_PVALB_CHC")

# Probability that eQTL is specific for each subclass
res = lapply(ctlist, function(query){

	include = grep(paste0("^", query), colnames(df_lfsr), value=TRUE)

	message(paste(include, collapse=", "))

	exclude = setdiff(colnames(df_lfsr), include)

	prob = dreamlet:::.compositePosteriorTest( 1 - df_lfsr, include = include, exclude = exclude, test="at least 1")

	res = data.frame(prob)
	colnames(res) = query
	res
})
res = bind_cols(res) %>%
		rownames_to_column("Gene") %>%
		as_tibble

ord.subclass = assay_order %>%
				grep("^subclass", ., value=TRUE) %>%
				gsub("^subclass / ", "", .) 

ids = colnames(res)[-1][!colnames(res)[-1] %in% ord.subclass]

for(query in ids){
	i = grep(paste0("^", query), ord.subclass)
 
 	ord.class <- append(ord.subclass, query, after=max(i))

 	cols[query] = mean_col(cols[ord.subclass[i]])
}

# Probability that eQTL is specific for each subclass
res_prob = res %>%
			pivot_longer(cols=!Gene, names_to = "subclass", values_to="prob") %>%
			mutate(ID = Gene) %>%
			mutate(Gene = gsub("^(\\S+)_.*$", "\\1", Gene)) 

res_prob %>%
	arrange(-prob) %>%
	# filter(prob > .1) %>%
	write_tsv(file="specific_eQTL_subclass.tsv.gz")

a = unique(res_prob$subclass)
ord.other = a[!a %in% ord.subclass]

fig = res_prob %>%
	group_by(subclass) %>%
	summarize( nSignif = sum(prob > 0.50)) %>%
	mutate(subclass = factor(subclass, c(ord.other, ord.subclass))) %>%
	ggplot(aes(nSignif, subclass, fill=subclass, label=nSignif)) +
		geom_bar(stat="identity") +
		theme_classic() +
		theme(aspect.ratio=1, legend.position = "none")  +
		scale_x_continuous(limits=c(0, 320), expand=c(0,0)) +
		scale_fill_manual(values = cols) +
		geom_text() +
		xlab("# genes with probability of specific eQTL > 50%") 

ggsave(fig, file="plots/specific_eQTL_subclass.pdf")



i = order(res_prob$prob, decreasing=TRUE)
head(res_prob[i,], 200) %>%
	data.frame

tail(res_prob[i,], 5000) %>%
	pull(Gene) %>%
	unique


# Gene heatmap
###############

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
		select(Gene, ID, subclass) %>%
		rename(CellType.specific = subclass)


# Get gene cell type pairs for 
# 1) Genes that coloc in any cell type and 
# 2) have a specific eQTL in any cell type
df_coloc %>%
	select(Gene, CellType) %>%
	rename(CellType.coloc = CellType) %>%
	distinct %>%
	inner_join(df1, relationship = "many-to-many") %>%
	filter(CellType.coloc != CellType.specific) %>%
	filter(!grepl("\\.\\*", CellType.specific)) %>%
	arrange(Gene) %>%
	write_tsv(file="merge_coloc_mashr.tsv")


df1 %>%
	filter(Gene == 'PFKFB2')

df_coloc %>%
	filter(Gene == 'PFKFB2')


df1 %>%
	filter(Gene == 'APP')

df_coloc %>%
	filter(Gene == 'APP')



# cell type specific eQTLs that also coloc with GWAS
res_spec_coloc = res_prob %>%
		mutate(CellType = factor(subclass, c(ord.other, ord.subclass))) %>%
		select(-subclass) %>%
		full_join(df_coloc, by=c("Gene"), relationship = "many-to-many") %>%
		filter(!is.na(Trait)) %>%
		filter(prob > 0.5) %>%
		arrange(Gene)


res_spec_coloc  %>% 
	filter(Gene == "APP")



res_prob %>% 
	filter(Gene == "APP") %>%
	arrange(-prob)



df_coloc %>% 
	filter(Gene == "APP")




# Reprt standard 
df_prob = df_lfsr %>%
	rownames_to_column("Gene") %>%
	pivot_longer(cols=!Gene, names_to = "subclass", values_to="lfsr") %>%
	mutate(CellType = factor(subclass,  c(ord.other, ord.subclass))) %>%
	select(-subclass)  %>%
	mutate(ID = Gene) %>%
	mutate(Gene = gsub("^(\\S+)_.*$", "\\1", Gene)) 

df_join = df_prob %>%
			left_join(df_coloc, relationship = "many-to-many") %>%
			# filter(!is.na(Trait)) %>%
			mutate(prob = 1 - lfsr)

df_overlap = df_join %>%
	filter(prob > 0.5) %>%
	filter(ppH4 > 0.5) %>%
	arrange( Gene, CellType ) %>%
	data.frame 

df_overlap %>%
	write_tsv("specific_and_coloc.tsv")



genes = c('NALCN', 'WNT5B', 'DOCK1', "GRIA1", 'APP', 'CACNA1C', 'PSD3', "INPP5D", "SERPINB1", "TLE4", res_spec_coloc$Gene)
genes = sort(unique(genes))

df_sub = df_join %>%	
	filter(Gene %in% genes) %>%
	mutate(Gene = factor(Gene, genes)) %>%
	mutate(CellType = factor(CellType, c(ord.other, ord.subclass))) %>%
	mutate(variant = gsub("^\\S*_(\\S+)$", "\\1", ID))


files3 = dir('/sc/arion/projects/CommonMind/hoffman/sceqtl/sumStats/', pattern="meta-.*-subclass.*parquet", full.names=TRUE)

dSet = open_dataset(files)

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
	# filter(Gene == "APP") %>% 
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



ord.class = assay_order %>%
				grep("^subclass", ., value=TRUE) %>%
				gsub("^subclass / ", "", .) 
		

fig = df_lfsr %>%
	rownames_to_column("Gene") %>%
	filter(Gene %in% genes) %>%
	pivot_longer(cols=!Gene, names_to = "class", values_to="prob") %>%
	mutate(class = factor(class, ord.class)) %>%
	mutate(Gene = factor(Gene, genes)) %>%
	ggplot(aes(Gene, class, fill = 1- prob)) +
		geom_tile() +
		theme_classic() + 
		coord_equal() +
		scale_fill_gradient(low="white", high="red", limits=c(0, 1)) +
		scale_x_discrete(guide = guide_axis(angle = 90)) 

ggsave(fig, file="plots/gene_heatmap.pdf")


cd /sc/arion/projects/CommonMind/hoffman/sceqtl/sumStats/
R 

suppressPackageStartupMessages({
library(tidyverse)
library(arrow)
library(ggplot2)
})


plotGenePanels = function(files, GENE, CT, SNP,window = 2e6, min.p=1e-7){

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
		ggplot(aes(Gene, name, fill=CPM)) +
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


# Get set of top variants
#-------------------------
# files1 = dir('/sc/arion/projects/CommonMind/hoffman/sceqtl/sumStats/', pattern="-bulk", full.names=TRUE)
# files2 = dir('/sc/arion/projects/CommonMind/hoffman/sceqtl/sumStats/', pattern="-class", full.names=TRUE)
files3 = dir('/sc/arion/projects/CommonMind/hoffman/sceqtl/sumStats/', pattern="meta-.*-subclass.*parquet", full.names=TRUE)

lst = list(c(GENE = "AUTS2",  CT = "Micro"),
		c(GENE = "CACNA1A", CT = "Micro"),
		c(GENE = "APOE", CT = "Micro"),
		c(GENE = "BDNF", CT = "EN_L2_3_IT"),
		c(GENE = "CLU", CT = "Micro"),
		c(GENE = "PDE4B", CT = "Oligo"),
		c(GENE = "SHANK2", CT = "Oligo"),
		c(GENE = "NXN", CT = "IN_PVALB_CHC"),
		c(GENE = "STAT3", CT = "IN_PVALB_CHC"),
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
	ggsave(res$fig.expr, file=file, height=10, width=4)
}
















