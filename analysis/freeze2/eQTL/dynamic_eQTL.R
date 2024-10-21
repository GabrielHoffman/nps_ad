


cd /sc/arion/projects/CommonMind/hoffman/sceqtl/sumStats/
R

library(ggplot2)
library(tidyverse)

assay_order = readRDS("assay_order.RDS")

ord = assay_order %>%
		grep("^class", ., value=TRUE) %>%
		gsub("^class / ", '', .)

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


files = system("ls /sc/arion/projects/CommonMind/zengb02/single_cell_eQTL_for_Gabriel_02212024/eQTL_detection_on_pearsonsubtype_level/dynamic_state_eQTL/Aging/*/analysis/merged_dynamic_eQTLs", intern=TRUE)

df_eqtl = lapply(files, function(file){

	ID = gsub("^.*Aging/(\\S+)/analysis.*$", "\\1", file )

	df = read_tsv(file, show_col_types=FALSE)
	df$CellType = ID
	df
	})
df_eqtl = bind_rows(df_eqtl)

df_eqtl %>%
	write_tsv("dynamic_eQTL.tsv")


files = system("ls /sc/arion/projects/CommonMind/zengb02/single_cell_eQTL_for_Gabriel_02212024/eQTL_detection_on_pearsonsubtype_level/dynamic_state_eQTL/Aging/*/analysis/significant_genes_dynamic_eQTL", intern=TRUE)

df_genes = lapply(files, function(file){

	ID = gsub("^.*Aging/(\\S+)/analysis.*$", "\\1", file )
	ID = ifelse(ID == 'Micro', "Immune", ID)

	df = read_tsv(file, col_names=FALSE, show_col_types=FALSE) %>%
			rename(Gene = X1)
	df$CellType = ID
	df
	})
df_genes = bind_rows(df_genes)

df_genes %>%
	write_tsv("dynamic_eQTL_significant_genes.tsv")


fig = df_genes %>%
		group_by(CellType) %>%
		summarize(nGenes = length(Gene)) %>%
		mutate(CellType = factor(CellType, ord)) %>%
		ggplot(aes(nGenes, CellType, fill=CellType, label=nGenes) ) +
			geom_bar(stat="identity") +
			theme_classic() +
			theme(aspect.ratio=1, legend.position="none") +
			scale_fill_manual(values=cols) +
			scale_x_continuous(name="# dynamic eQTLs", limits=c(0, 1500), expand=c(0,0)) +
			geom_text(color="black")

ggsave(fig, file="plots/dynamic_eQTL_counts.pdf")


# Coloc
#######

cd /Users/gabrielhoffman/Library/CloudStorage/Dropbox/projects/PsychAD3.0/supplementary_figures/dynamic_eQTL_colocalization
R

library(ggplot2)
library(tidyverse)

file = "colocalization_results_for_dynamic_eQTL.tsv"
df = read_tsv(file) %>%
	mutate(gene = factor(gene, rev(sort(unique(gene))))) 
colnames(df) = c("CellType", "Gene", "Trait", "ppH4")	

fig = df %>%
	ggplot( aes(Trait, Gene, color=ppH4, size=ppH4, label=ifelse(ppH4 > 0.8, "x", ''))) +
	geom_point() +
	theme_classic() +
	# coord_equal() +
	theme(plot.title = element_text(hjust = 0.5)) +
	scale_size_area() +
	scale_color_gradient(low="white", high="red", limits=c(0, 1)) +
	geom_text(color="black", vjust=0.5, hjust=0.5) +
	facet_wrap(~CellType, nrow=1, scales="free") 
ggsave(fig, file="dynamic_coloc_all_v1.pdf", height=4, width=12)



fig = df %>%
	filter(ppH4 > .5) %>%
	ggplot(aes(CellType, Gene, color=ppH4, size=ppH4, label=ifelse(ppH4 > 0.8, "x", ''))) +
		geom_point() +
		theme_classic() +
		# coord_equal() +
		theme(plot.title = element_text(hjust = 0.5)) +
		facet_wrap(~Trait, nrow=1, scales="free") +
		scale_size_area() +
		scale_color_gradient(low="white", high="red", limits=c(0, 1)) +
		geom_text(color="black", vjust=0.5, hjust=0.5)
ggsave(fig, file="dynamic_coloc_all_v2.pdf", height=4, width=9)



df %>%	
	filter(ppH4 > .01) %>%
	mutate(ylabel = paste(Trait, '/', Gene)) %>%
	ggplot(aes(CellType, ylabel, color=ppH4, size=ppH4, label=ifelse(ppH4 > 0.8, "x", ''))) +
		geom_point() +
		theme_classic() +
		coord_equal() +
		theme(plot.title = element_text(hjust = 0.5)) +
		scale_size_area() +
		scale_color_gradient(low="white", high="red", limits=c(0, 1)) +
		geom_text(color="black", vjust=0.5, hjust=0.5)

ord = rev(c("Micro", "Astro", "OPC", "Oligo", "EN", "IN"))
cols = c( '#1A9993','#197EC0','#D5E4A2','#FED439','#D2AF81','#C80813')
names(cols) = ord

# count
fig = df %>%
	filter(ppH4 > .5)  %>%
		group_by(CellType, Trait) %>%
		summarize(nGenes = length(Gene)) %>%
		mutate(CellType = factor(CellType, ord))  %>%
		ggplot(aes(CellType, nGenes, fill=CellType, label=nGenes) ) +
			geom_bar(stat="identity", position = "dodge2") +
			theme_classic() +
			theme(aspect.ratio=1, legend.position="none") +
			scale_fill_manual(values=cols) +
			scale_y_continuous(name="# dynamic eQTLs", limits=c(0, 8.3), expand=c(0,0)) +
			geom_text(color="black") +
			facet_wrap(~Trait, ncol=1)

file = "coloc_dynamic_total.pdf"
ggsave(fig, file=file)


expression ~ trajectory + SNP + trajectory*SNP

# coloc with _only_ dynamic eQTL
################################

file = "colocalization_results_for_dynamic_but_not_standard_eQTL.tsv"
df = read_tsv(file, col_names=c("CellType", "Gene", "Trait", "ppH4")) %>% 
	mutate(Gene = factor(Gene, rev(sort(unique(Gene)))))


fig = df %>%
	mutate(ylabel = paste(Trait, '/', Gene)) %>%
	mutate(ylabel = factor(ylabel, rev(unique(sort(ylabel))))) %>%
	ggplot(aes(CellType, ylabel, color=ppH4, size=ppH4, label=ifelse(ppH4 > 0.8, "x", ''))) +
		geom_point() +
		theme_classic() +
		theme(plot.title = element_text(hjust = 0.5)) +
		scale_size_area() +
		scale_color_gradient(low="white", high="red", limits=c(0, 1)) +
		geom_text(color="black", vjust=0.5, hjust=0.5)

ggsave(fig, file="dynamic_coloc_highlight.pdf", height=5, width=3.3)











