


cd /sc/arion/projects/CommonMind/hoffman/sceqtl/sumStats/
R

library(ggplot2)
library(tidyverse)

files = system("ls /sc/arion/projects/CommonMind/zengb02/single_cell_eQTL_for_Gabriel_02212024/eQTL_detection_on_pearsonsubtype_level/dynamic_state_eQTL/Aging/*/analysis/merged_dynamic_eQTLs", intern=TRUE)

df_eqtl = lapply(files, function(file){

	ID = gsub("^.*Aging/(\\S+)/analysis.*$", "\\1", file )

	df = read_tsv(file)
	df$CellType = ID
	df
	})
df_eqtl = bind_rows(df_eqtl)


# Coloc
#######

cd /Users/gabrielhoffman/Library/CloudStorage/Dropbox/projects/PsychAD3.0/supplementary_figures/dynamic_eQTL_colocalization


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




# coloc with _only_ dynamic eQTL
################################

file = "colocalization_results_for_dynamic_but_not_standard_eQTL.tsv"
df = read_tsv(file, col_names=c("CellType", "Gene", "Trait", "ppH4")) %>% 
	mutate(Gene = factor(Gene, rev(sort(unique(Gene)))))


fig = df %>%
	mutate(ylabel = paste(Trait, '/', Gene)) %>%
	ggplot(aes(CellType, ylabel, color=ppH4, size=ppH4, label=ifelse(ppH4 > 0.8, "x", ''))) +
		geom_point() +
		theme_classic() +
		theme(plot.title = element_text(hjust = 0.5)) +
		scale_size_area() +
		scale_color_gradient(low="white", high="red", limits=c(0, 1)) +
		geom_text(color="black", vjust=0.5, hjust=0.5)

ggsave(fig, file="dynamic_coloc_all_v2.pdf", height=4, width=9)









