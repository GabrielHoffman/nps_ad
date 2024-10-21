
# trans and cis heritability
# /sc/arion/projects/CommonMind/zengb02/single_cell_eQTL_for_Gabriel_02212024/eQTL_detection_on_pearson*_level/heritability_analysis/estimate_heritability_explained/merged_expression_heritability_*_level

# Write eQTL summary stats to PARQUET
######################################

library(parallel)
library(tidyverse)
library(arrow)

# MSSM
# files = system("ls /sc/arion/projects/CommonMind/zengb02/single_cell_eQTL_for_Gabriel_02212024/eQTL_detection_on_pearson*_level/cis-eQTL_detection/eQTL_results/eQTL_result_MSSM_*_chr*.gz", intern=TRUE)

# meta analysis eQTLs
files = system("ls /sc/arion/projects/CommonMind/zengb02/single_cell_eQTL_for_Gabriel_02212024/eQTL_detection_on_pearson*_level/cis-eQTL_detection/eQTL_results/meta_results/eQTL_result_merged_*chr*.gz", intern=TRUE)

ptrn = "^eQTL_result_([a-zA-Z]+)_((\\S+)+)_(\\S+).gz$"
info = lapply(files, function(file){
	Source = gsub(ptrn, "\\1", basename(file))
	Source = ifelse(Source == "merged", "meta", Source)
	Level = gsub(".*pearson([a-zA-Z]+)_level.*", "\\1", file)
	CellType = gsub(ptrn, "\\2", basename(file))
	Chrom = gsub(ptrn, "\\3", basename(file))

	tibble(file, Source, Level, CellType, Chrom)
})
info = bind_rows(info)
info$key = with(info, paste(Source, CellType, Level, sep='-'))

# table(info$Source)
# table(info$CellType, info$Level)
# table(info$Chrom)

# get eGenes
#---------

# MSSM significant genes
# files = system("ls /sc/arion/projects/CommonMind/zengb02/single_cell_eQTL_for_Gabriel_02212024/eQTL_detection_on_pearson*_level/cis-eQTL_detection/analysis/significant_eQTL*/significant_peaks_MSSM_*", intern=TRUE)

# df_QTL = mclapply(files, function(file){
# 	ptrn = "^significant_peaks_([a-zA-Z]+)_((\\S+)+)$"
# 	Source = gsub(ptrn, "\\1", basename(file))
# 	Level = gsub(".*pearson([a-zA-Z]+)_level.*", "\\1", file)
# 	CellType = gsub(ptrn, "\\2", basename(file))
	
# 	df = read_tsv(file, col_names=cn, progress=FALSE, show_col_types = FALSE) %>%
# 		select(gene, top_variant, z_score, FDR_gene)

# 	df$Source = Source
# 	df$Level = Level
# 	df$CellType = CellType
# 	df
# }, mc.cores=24)
# df_QTL = bind_rows(df_QTL)


# meta significant genes
files = system("ls /sc/arion/projects/CommonMind/zengb02/single_cell_eQTL_for_Gabriel_02212024/eQTL_detection_on_pearson*_level/cis-eQTL_detection/meta_eQTL_detection/adjust_combined_result_*_allChr", intern=TRUE)

cn = c("gene", "top_variant", "z_score", "adjust_p_gene_level1", "adjust_p_gene_level2", "adjust_p_gene_level3", 'FDR_gene')

df_QTL = mclapply(files, function(file){ 
	Source = "meta"
	Level = gsub(".*pearson([a-zA-Z]+)_level.*", "\\1", file)
	CellType = gsub("^adjust_combined_result_(\\S+)_allChr$", "\\1", basename(file))
	
	df = read_tsv(file, col_names=cn, progress=FALSE, show_col_types = FALSE) %>%
		select(gene, top_variant, z_score, FDR_gene)

	df$Source = Source
	df$Level = Level
	df$CellType = CellType
	df
}, mc.cores=24)
df_QTL = bind_rows(df_QTL)

write_parquet(df_QTL, "df_QTL.parquet")



# setup to read files
cn = c('chr', 'start', 'end', 'variant', 'nonassessed_allele', 'assessed_allele', 'gene', 'beta', 'se', 'z')

# for each 
res = mclapply( unique(info$key), function(keyValue){
	message(keyValue)
	# get file paths
	idx = which(info$key == keyValue)
	files = info$file[idx]
	files = files[file.size(files) > 100]

	# read files
	df = read_tsv(files, col_names=cn, progress=FALSE, show_col_types = FALSE, num_threads=1)

	if(nrow(df) == 0) return(NA)

	df$Source = info$Source[idx][1]
	df$Level = info$Level[idx][1]
	df$CellType = info$CellType[idx][1]

	outFile = paste0('/sc/arion/projects/CommonMind/hoffman/sceqtl/sumStats/', keyValue, ".parquet")

	res = df %>% 
		select(Source, Level, CellType, all_of(cn)) %>%
		arrange(gene, -abs(z)) %>%
		group_by(gene) %>%
		mutate(isTopHit = (variant == first(variant))) %>%
		write_parquet(outFile)
	NA
}, mc.cores=24, mc.preschedule=FALSE)

# Read data for analysis
########################

cd /sc/arion/projects/CommonMind/hoffman/sceqtl/sumStats/
R 

suppressPackageStartupMessages({
library(tidyverse)
library(arrow)
library(ggplot2)
library(parallel)
library(cowplot)
library(ggrepel)
})

# genome-wide significant eGenes
file = "/sc/arion/projects/CommonMind/hoffman/sceqtl/sumStats/df_QTL.parquet"
df_QTL = read_parquet( file ) %>%
			filter( !(Level == 'subclass' & CellType %in% c("IN", "EN", "Immune",  "Mural")))

# Get set of top variants
#-------------------------
files1 = dir('/sc/arion/projects/CommonMind/hoffman/sceqtl/sumStats/', pattern="-bulk", full.names=TRUE)
files2 = dir('/sc/arion/projects/CommonMind/hoffman/sceqtl/sumStats/', pattern="-class", full.names=TRUE)
files3 = dir('/sc/arion/projects/CommonMind/hoffman/sceqtl/sumStats/', pattern="-subclass", full.names=TRUE)
files4 = dir('/sc/arion/projects/CommonMind/hoffman/sceqtl/sumStats/', pattern="-subtype", full.names=TRUE)
files = c(files1, files2, files3, files4)
dSet = open_dataset(files)

topVariants = dSet %>%
				filter(isTopHit) %>%
				select(variant) %>%
				distinct %>%
				collect
# saveRDS(topVariants, file="topVariants.RDS")
topVariants = readRDS("topVariants.RDS")


# Extract eQTL results from these top variants in any cell type
df_all = dSet %>%
			select(Source, Level, CellType, variant, gene, beta, se, isTopHit) %>%
			inner_join(topVariants, by="variant") %>%
			inner_join(df_QTL %>% select(-top_variant, -z_score), by=c("gene", 'Source', 'Level', 'CellType'))  %>%
			mutate(Dataset = paste0(Level, ' / ', CellType)) %>%
			collect
# saveRDS(df_all, file="df_all.RDS")
df_all = readRDS("df_all.RDS")


# Cell frequency
library(dreamlet)
library(parallel)
library(tidyverse)
files = c(class = "/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/pseudobulk/MSSM_2024-02-01_16_17_PB_SubID_class.RDS", 
		subclass = "/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/pseudobulk/MSSM_2024-02-01_16_17_PB_SubID_subclass.RDS", 
		subtype = "/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/pseudobulk/MSSM_2024-02-01_16_17_PB_SubID_subtype.RDS")
df_frac = mclapply(names(files), function(id){

	library(dreamlet)

	pb = readRDS(files[id])

	count = colSums(cellCounts(pb))
	tibble(Dataset = paste0(id, " / ", names(count)), 
		Fraction = count / sum(count))
	}, mc.cores=3)
df_frac = bind_rows(df_frac)
# saveRDS(df_frac, file="df_frac.RDS")
df_frac = readRDS("df_frac.RDS")

# get colors
file = "/sc/arion/projects/psychAD/NPS-AD/freeze2_proc/230816_PsychAD_capstone_F1/230921_PsychAD_color_palette.csv"
df_colors = read_csv(file) %>%
				select(category, name, color_hex) %>%
				filter(category %in% c("class", "subclass", "subtype")) %>%
				mutate(Dataset = paste0(category, " / ", name)) %>%
				select(Dataset, color = color_hex) %>%
				bind_rows(tibble(Dataset = "bulk / bulk", color = "grey"))
cols = df_colors$color
names(cols) = df_colors$Dataset

# Get read counts per subject
df_reads = readRDS("df_counts.RDS") %>%
	mutate(Dataset = paste0(Level,' / ', celltype))

# Cell type hierarchy
#####################
library(ape)
library(dendextend)
library(ggplot2)
library(ggtree)

plot_tree_simple = function(tree, xmax.scale=1.5){

    fig.tree = ggtree(tree, branch.length = "none", ladderize=FALSE) + 
               geom_tiplab(color = "black", size=4, hjust=0, offset=.2) +
               theme(legend.position="top left", plot.title = element_text(hjust = 0.5))

    # get default max value of x-axis
    xmax = layer_scales(fig.tree)$x$range$range[2]

    # increase x-axis width
    fig.tree + xlim(0, xmax*xmax.scale) 
}

# load
files = system("ls /sc/arion/projects/psychAD/NPS-AD/freeze2_proc/231109_PsychAD_capstone_F2/tree_*_um.nwk", inter=TRUE)

plotList = lapply( files, function(file){
	tree = read.tree(file)

	tree = drop.tip(tree, "Adaptive")
	tree = drop.tip(tree, "EN_L5_ET")

	plot_tree_simple(as.phylo(tree), xmax.scale=1.5) + theme(legend.position="bottom")
})

pdf("plots/trees.pdf")
plotList
dev.off()


tree = read.tree(files[1])
tree = drop.tip(tree, "Adaptive")
tree = drop.tip(tree, "EN_L5_ET")
tree$tip.label[] = ''
fig = plot_tree_simple(as.phylo(tree), xmax.scale=1.5) + theme(legend.position="bottom")
ggsave(fig, file="plots/trees_1.pdf", height=.5*4, width=.17*4)


tree = read.tree(files[2])
tree = drop.tip(tree, "Adaptive")
tree = drop.tip(tree, "EN_L5_ET")
tree$tip.label[] = ''
fig = plot_tree_simple(as.phylo(tree), xmax.scale=1.5) + theme(legend.position="bottom")
ggsave(fig, file="plots/trees_2.pdf", height=1.76*4, width=.17*4)



assay_order = lapply(plotList, ggtree::get_taxa_name)

assay_order = c("bulk / bulk", paste0("class / ", assay_order[[1]]),
	paste0("subclass / ", assay_order[[2]]),
	paste0("subtype / ", assay_order[[3]]))
# saveRDS(assay_order, file="assay_order.RDS")



# plot of # eGenes
#-----------------

# Upset: class
res_gene = df_QTL %>%
			filter(Level == "class") 

qtlLst = lapply(unique(res_gene$CellType), function(CT){
	res_gene %>%
		filter(FDR_gene < 0.05) %>%
		filter(CellType == CT) %>%
		pull(gene) %>%
		unique
	})
names(qtlLst) = unique(res_gene$CellType)

assay_order = readRDS("assay_order.RDS")
ord = assay_order %>%
		grep("^class", ., value=TRUE) %>%
		gsub("^class / ", '', .)

pdf("plots/standard_upset_class.pdf", height=5, width=7)
upset(fromList(qtlLst[ord]), order.by = "freq", nsets=8, keep.order=TRUE, sets=ord, nintersects=15)
dev.off()

# Upset: subclass
res_gene = df_QTL %>%
			filter(Level == "subclass") 

qtlLst = lapply(unique(res_gene$CellType), function(CT){
	res_gene %>%
		filter(FDR_gene < 0.05) %>%
		filter(CellType == CT) %>%
		pull(gene) %>%
		unique
	})
names(qtlLst) = unique(res_gene$CellType)

assay_order = readRDS("assay_order.RDS")
ord = assay_order %>%
		grep("^subclass", ., value=TRUE) %>%
		gsub("^subclass / ", '', .)

ord = ord[ord %in% names(qtlLst)]

pdf("plots/standard_upset_subclass.pdf", height=9, width=7)
upset(fromList(qtlLst[ord]), order.by = "freq", nsets=8, keep.order=TRUE, sets=ord, nintersects=15)
dev.off()






# For Oligo, keep from top Level
df4 = df_QTL %>%
	filter(FDR_gene < 0.05) %>%
	group_by(Source, Level, CellType) %>%
	summarize(nGenes = length(gene))  %>%
	mutate(Dataset = paste0(Level,' / ', CellType)) %>%
	mutate(isNeuron = grepl("EN|IN", Dataset)) %>%
	mutate(Type = gsub("^(.*) /.*$", "\\1", Dataset)) %>%
	mutate(Celltype = gsub("^(.*) / (.*)$", "\\2", Dataset)) %>%
	left_join(df_frac, by="Dataset") %>%
	mutate(Level = factor(Level, c("bulk", "class", "subclass", "subtype"))) %>%
	left_join(df_reads, by="Dataset") %>%
	mutate(Level = Level.x) %>%
	arrange(Level, CellType)  


plot_eGene_count = function( df4, lvl){
	# get Rsq value
	f = function(data){
		fit = lm(nGenes ~ Fraction, data = data)
		rsq = summary(fit)$r.sq
		data.frame(rsq = rsq, p = coef(summary(fit))[2,4])
	}
	df_rsq = df4 %>%
		filter(Level %in% lvl) %>%
		# group_by(isNeuron) %>%
		group_modify(~f(.))

	# i = which(df_rsq$isNeuron)
	i = 1
	df_rsq$text = NA
	df_rsq$text[i] = paste("Rsq:", format(df_rsq$rsq[i], digits=3), "\np:", format(df_rsq$p[i], digits=3))
	# df_rsq$text[-i] = paste("Non-neuron:", format(df_rsq$rsq[-i], digits=3))

	# xlimit = range(100*df4$Fraction, na.rm=TRUE)
	ymax = df4 %>%
		filter(Level %in% lvl) %>%
		pull(nGenes) %>%
		max

	fig = df4 %>%
		filter(Level %in% lvl) %>%
		filter(! duplicated(Celltype)) %>%	
		ggplot(aes(100*Fraction, nGenes, color=Dataset, group=1, shape=Type)) + #isNeuron
			geom_point(size=3) +
			theme_classic() +
			# facet_wrap(~Level) +
			scale_color_manual(name = "Cell type", values=cols) +
			theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5), legend.position="none") +
			# scale_x_log10(limit=xlimit) +
			scale_x_log10() +
			geom_smooth(method = "lm", se=FALSE) +
			scale_y_continuous(limits=c(0,ymax*1.04), expand=c(0,0)) +
			xlab("Cell type abundance (%)") +
			ylab("# of eGenes")  +
			geom_text_repel(aes(label=Celltype), box.padding=.4, max.overlaps=10, size=4) +
			annotate("text", x= 1, y = 5000, label=df_rsq$text) +
			ggtitle(paste(lvl, collapse=', '))
}

df5 = df4 %>%
		filter(Type != "bulk")

fig1 = plot_eGene_count( df5, "class")
fig2 = plot_eGene_count( df5, "subclass")
fig3 = plot_eGene_count( df5, "subtype")
fig = plot_grid(fig1, fig2, fig3, nrow=1)
ggsave(fig, file="plots/eGenes_v1.pdf", height=5, width=12)


fig1 = plot_eGene_count( df4, c("class", "subclass"))
fig2 = plot_eGene_count( df4, "subtype")
fig = plot_grid(fig1, fig2, nrow=1)
# ggsave(fig, file="~/www/test.pdf", height=5, width=8)
ggsave(fig, file="plots/eGenes_v2.pdf", height=5, width=8)

# eGenes vs total read count
plot_eGene_count = function( df4, lvl){
	# get Rsq value
	f = function(data){
		fit = lm(nGenes ~ meanReadCount, data = data)
		rsq = summary(fit)$r.sq
		data.frame(rsq = rsq, p = coef(summary(fit))[2,4])
	}
	df_rsq = df4 %>%
		filter(Level %in% lvl) %>%
		group_modify(~f(.))

	df_rsq$text = NA
	df_rsq$text = paste("Rsq:", format(df_rsq$rsq, digits=3), "\np:", df_rsq$p)

	ymax = df4 %>%
		filter(Level %in% lvl) %>%
		pull(nGenes) %>%
		max

	fig = df4 %>%
		filter(Level %in% lvl) %>%
		filter(! duplicated(Celltype)) %>%	
		ggplot(aes(meanReadCount, nGenes, color=Dataset, shape=Type, group=1)) +
			geom_point(size=3) +
			theme_classic() +
			# facet_wrap(~Level) +
			scale_color_manual(name = "Cell type", values=cols) +
			theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5), legend.position="none") +
			geom_smooth(method = "lm", se=FALSE) +
			scale_y_continuous(limits=c(0,ymax*1.04), expand=c(0,0)) +
			scale_x_log10() +
			xlab("Mean reads per Subject") +
			ylab("# of eGenes")  +
			geom_text_repel(aes(label=Celltype), box.padding=.4, max.overlaps=10, size=4) +
			annotate("text", x= 2000000, y = 2000, label=df_rsq$text) +
			ggtitle(paste(lvl, collapse=', '))
}

fig1 = plot_eGene_count( df4, c("class"))
fig2 = plot_eGene_count( df4, "subclass")
fig3 = plot_eGene_count( df4, "subtype")
fig = plot_grid(fig1, fig2, fig3, nrow=1)
ggsave(fig, file="plots/eGenes_readCount_v1.pdf", height=5, width=12)


fig1 = plot_eGene_count( df4, c("class", "subclass"))
fig2 = plot_eGene_count( df4, "subtype")
fig = plot_grid(fig1, fig2, nrow=1)
# ggsave(fig, file="~/www/test.pdf", height=5, width=8)
ggsave(fig, file="plots/eGenes_readCount_v2.pdf", height=5, width=8)


# Barplot of # eGenes
# ymax = max(df4$nGenes) * 1.02
ymax = 11000
fig = df4 %>%
	mutate(Dataset = factor(Dataset, assay_order)) %>%
	# filter(is.na(Fraction))
	filter(Type %in% c("class", "subclass")) %>%
	ggplot(aes(Dataset, nGenes, fill=Dataset, label=nGenes)) +
		geom_bar(stat="identity") +
		theme_classic() +
		coord_flip() +
		scale_fill_manual(name = "Cell type", values=cols) +
		theme(aspect.ratio=3, legend.position="none") +
		scale_y_continuous(limits=c(0, ymax), expand=c(0,0)) +
		geom_text()
ggsave(fig, file="plots/eGenes_v3.pdf", height=12, width=8)

# scp sklar1:"/sc/arion/projects/CommonMind/hoffman/sceqtl/sumStats/plots/*" .

# Correlation with bulk eQTLs
#----------------------------

df_all = readRDS("df_all.RDS")

CT1 = "bulk / bulk"

a = df_all %>%
	filter(Dataset == CT1)

ids = unique(df_all$Dataset)
ids = ids[ids != CT1]

df2 = lapply(ids, function(CT2){
	message(CT2)

	b = df_all %>%
		filter(Dataset == CT2)

	df_join1 = inner_join(a,b, 
					by = c( "Source", "variant", "gene")) %>% 
					filter(isTopHit.x, FDR_gene.x < 0.05) %>%
					mutate(discovery = CT1) 

	df_join2 = inner_join(a,b, 
					by = c( "Source", "variant", "gene")) %>% 
					filter(isTopHit.y, FDR_gene.y < 0.05) %>%
					mutate(discovery = CT2)

	df = bind_rows(df_join1, df_join2) %>%
			select(gene, variant, Dataset.x, beta.x, FDR_gene.x, Dataset.y,  beta.y, FDR_gene.y, discovery)   %>%
		mutate(agree = sign(beta.x) == sign(beta.y))  %>%
		mutate(c.x = beta.x*sign(beta.y),
				c.y = beta.y*sign(beta.x)) %>%
		mutate(discovery = factor(discovery, c(CT1, CT2)))

	df %>%
		group_by(discovery) %>%
		summarise(cor = cor(beta.x, beta.y)) %>%
		mutate(Dataset = CT2)
})
df2 = bind_rows(df2)

df2$isNeuron = grepl("EN|IN", df2$Dataset)
df2$Type = gsub("^(.*) /.*$", "\\1", df2$Dataset)
df2$Celltype = gsub("^(.*) / (.*)$", "\\2", df2$Dataset)


df3 = df2 %>%
	arrange(cor) %>%
	mutate(discovery = (discovery != 'bulk / bulk'))  %>%
	inner_join(df_frac, by='Dataset') %>%
	filter(!discovery) %>%
	filter(! Dataset %in% c("subclass / Astro", "subclass / Oligo", "subclass / Endo", "subclass / OPC")) 

# get Rsq value
f = function(data){
	fit = lm(cor ~ Fraction, data = data)
	rsq = summary(fit)$r.sq
	data.frame(rsq = rsq, p = coef(summary(fit))[2,4])
}
df_rsq = df3 %>%
	group_by(isNeuron) %>%
	group_modify(~f(.))

i = which(df_rsq$isNeuron)
df_rsq$text = NA
df_rsq$text[i] = paste("Neuron:", format(df_rsq$rsq[i], digits=3), "\n")
df_rsq$text[-i] = paste("Non-neuron:", format(df_rsq$rsq[-i], digits=3))


fig = df3 %>%
	ggplot(aes(100*Fraction, cor, group=isNeuron, color=Dataset, shape=Type)) +
		geom_smooth(method = "lm", se=FALSE) +	
		geom_point(size=3) +	
		theme_classic() +
		theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5), legend.position="none") +
		ylab("Correlation of estimated effect size comapred to bulk")  +
		xlab("Cell type abundance (%)") +
		scale_x_log10() +
		scale_y_continuous(limits=c(0,1), expand=c(0,0)) +
		scale_color_manual(name = "Cell type", values=cols) +
		geom_text_repel(aes(label=Celltype), box.padding=.4, max.overlaps=15, size=4) +
		annotate("text", x= 10,y = .25, label=df_rsq$text)

# ggsave(fig, file="~/www/test.pdf", height=5, width=5)
ggsave(fig, file="plots/eGenes_scatterplot.pdf", height=5, width=5)



# Effect size
##############


df_bulk = df_all %>%
			filter(Dataset == "bulk / bulk", isTopHit)

# all variant/genes that are significant in Bulk
df2 = df_all %>%
		inner_join( df_bulk, by = c("variant", "gene"))

fig = df2 %>%
	ggplot(aes(Dataset.x, beta.x)) +
		geom_violin() +
		coord_flip()
ggsave(fig, file="~/www/test.png")

# all variant/genes that are significant in cell type
df_ct = df_all %>%
			filter(isTopHit) %>%
			inner_join(df_all %>% filter(Dataset == "bulk / bulk"), by = c("variant", "gene"))

fig = df_ct %>%
	ggplot(aes(Dataset.x, abs(beta.x))) +
		geom_violin() +
		geom_boxplot(width=.1) +
		coord_flip()
ggsave(fig, file="~/www/test.png")


fig = df_ct %>%
	ggplot(aes(Dataset.x, beta.x - beta.y)) +
		geom_violin() +
		geom_boxplot(width=.1) +
		coord_flip()
ggsave(fig, file="~/www/test.png")


# Microglia vs bulk
df2 = df_all %>%
			filter(Dataset == "class / Immune", isTopHit)%>%
			inner_join(df_all %>% filter(Dataset == "bulk / bulk"), by = c("variant", "gene"))

fig = df2 %>%
	# ggplot(aes(beta.x, beta.y)) +
	ggplot(aes(abs(beta.x), abs(beta.y))) +
		geom_point() +
		theme_bw() +
		geom_abline(intercept=0, slope=1) +
		xlab("Immune") +
		ylab("bulk")
ggsave(fig, file="~/www/test.png")


# effect size for all cell types
fig = df_all %>%
		filter(isTopHit) %>%
		ggplot(aes(abs(beta), Dataset)) +
			geom_violin() +
			geom_boxplot(width=.1) +
			theme_bw() 
ggsave(fig, file="~/www/test.png")






df2 = df_all %>%
			filter(variant %in% topVariants) 



df3 = df2 %>%
		inner_join(df2 %>% filter(Dataset == "bulk / bulk"), by = c("variant", "gene"))



fig = df3 %>%
	ggplot(aes(Dataset.x, beta.x - beta.y)) +
		geom_violin() +
		geom_boxplot(width=.1) +
		coord_flip()
ggsave(fig, file="~/www/test.png")

topVariants



df_all %>%
	group_by(Dataset) %>%
	summarize(n = length(unique(.$gene)))




df_all %>%
	filter(isTopHit) %>%
	group_by(Dataset) %>%
	summarize(n = n())

df_all %>%
	filter(isTopHit) %>%
	ggplot(aes(Dataset, abs(beta))) +
		geom_violin()






# OLD
##########


df2 %>%
	arrange(cor) %>%
	mutate(discovery = (discovery != 'bulk'))  %>%
	inner_join(df_frac, by='Dataset') %>%
	ggplot(aes(Dataset, cor, fill=Dataset)) +
		geom_bar(stat="identity") +
		theme_classic() +
		theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5), legend.position="none") +
		coord_flip() +
		scale_y_continuous(limits=c(0,1), expand=c(0,0)) +
		facet_wrap(~ discovery) +
		ggtitle("Discovery in Bulk") +
		ylab("Correlation in estimated effect size") 










CT1 = "bulk"
CT2 = "OPC"

a = df_all %>%
	filter(CellType == CT1)

b = df_all %>%
	filter(CellType == CT2)

df_join = inner_join(a,b, by=c( "Source", "variant", "gene")) 				%>% filter(isTopHit.x | isTopHit.y)

fig1 = df_join %>%
	filter(isTopHit.x, FDR_gene.x < 0.05) %>%
	mutate(agree = sign(beta.x) == sign(beta.y))  %>%
	ggplot(aes(beta.x, beta.y, color=agree)) +
		geom_point() +
		theme_classic() +
		geom_abline(color="red") +
		xlab(CT1) +
		ylab(CT2)+
		coord_equal()


df_join = inner_join(a,b, by=c( "Source", "variant", "gene")) 				%>% filter(isTopHit.x | isTopHit.y)

fig2 = df_join %>%
	filter(isTopHit.y, FDR_gene.y < 0.05) %>%
	mutate(agree = sign(beta.x) == sign(beta.y))  %>%
	ggplot(aes(beta.x, beta.y, color=agree)) +
		geom_point() +
		theme_classic() +
		geom_abline(color="red") +
		xlab(CT1) +
		ylab(CT2) +
		coord_equal()

plot_grid(fig1, fig2)


df = df_join %>%
	filter(isTopHit.x, FDR_gene.x < 0.05) %>%
	mutate(agree = sign(beta.x) == sign(beta.y))  %>%
	mutate(c.x = beta.x*sign(beta.y),
			c.y = beta.y*sign(beta.x)) 

fig1 = df %>% 
	ggplot(aes(c.x,c.y, color=agree)) +
		geom_point() +
		theme_classic() +
		geom_abline(color="red") +
		xlab(CT1) +
		ylab(CT2)+
		coord_equal()


df = df_join %>%
	filter(isTopHit.y, FDR_gene.y < 0.05) %>%
	mutate(agree = sign(beta.x) == sign(beta.y))  %>%
	mutate(c.x = beta.x*sign(beta.y),
			c.y = beta.y*sign(beta.x)) 

fig2 = df %>% 
	ggplot(aes(c.x,c.y, color=agree)) +
		geom_point() +
		theme_classic() +
		geom_abline(color="red") +
		xlab(CT1) +
		ylab(CT2)+
		coord_equal()

plot_grid(fig1, fig2)





CT1 = "bulk"
CT2 = "Immune"

a = df_all %>%
	filter(CellType == CT1)

b = df_all %>%
	filter(CellType == CT2)


df_join1 = inner_join(a,b, 
				by = c( "Source", "variant", "gene")) %>% 
				filter(isTopHit.x, FDR_gene.x < 0.05) %>%
				mutate(discovery = CT1) 

df_join2 = inner_join(a,b, 
				by = c( "Source", "variant", "gene")) %>% 
				filter(isTopHit.y, FDR_gene.y < 0.05) %>%
				mutate(discovery = CT2)

df = bind_rows(df_join1, df_join2) %>%
		select(gene, variant, CellType.x, beta.x, FDR_gene.x, CellType.y,  beta.y, FDR_gene.y, discovery)   %>%
	mutate(agree = sign(beta.x) == sign(beta.y))  %>%
	mutate(c.x = beta.x*sign(beta.y),
			c.y = beta.y*sign(beta.x)) 


df %>%
	ggplot(aes(discovery, c.y - c.x , fill=discovery)) +
		geom_violin() +
		# geom_boxplot(width=.1) +
		geom_hline(yintercept=0) +
		ggtitle(paste(CT2, "-", CT1)) +
		theme_classic() +
		theme(legend.position="none") +
		ylim(-1,1) +
		facet_wrap(~agree)






CT1 = "bulk"

a = df_all %>%
	filter(CellType == CT1)

ids = unique(df_all$CellType)
ids = ids[ids != "bulk"]

figList = lapply(ids, function(CT2){
	message(CT2)

	b = df_all %>%
		filter(CellType == CT2)

	df_join1 = inner_join(a,b, 
					by = c( "Source", "variant", "gene")) %>% 
					filter(isTopHit.x, FDR_gene.x < 0.05) %>%
					mutate(discovery = CT1) 

	df_join2 = inner_join(a,b, 
					by = c( "Source", "variant", "gene")) %>% 
					filter(isTopHit.y, FDR_gene.y < 0.05) %>%
					mutate(discovery = CT2)

	df = bind_rows(df_join1, df_join2) %>%
			select(gene, variant, CellType.x, beta.x, FDR_gene.x, CellType.y,  beta.y, FDR_gene.y, discovery)   %>%
		mutate(agree = sign(beta.x) == sign(beta.y))  %>%
		mutate(c.x = beta.x*sign(beta.y),
				c.y = beta.y*sign(beta.x)) %>%
		mutate(discovery = factor(discovery, c(CT1, CT2)))

	df %>%
		ggplot(aes(discovery, c.y - c.x , fill=discovery)) +
			geom_violin() +
			# geom_boxplot(width=.1) +
			geom_hline(yintercept=0) +
			ggtitle(paste(CT2, "-", CT1)) +
			theme_classic() +
			theme(legend.position="none", aspect.ratio=1) +
			ylim(-1,1)  +
		facet_wrap(~agree)
})

fig = plot_grid(plotlist = figList, ncol=2)

ggsave(fig, file="~/www/test.png", height=10, width=9)








fig = plot_grid(plotlist = figList, ncol=2)

ggsave(fig, file="~/www/test.png", height=10, width=9)







