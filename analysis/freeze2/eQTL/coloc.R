



cd /sc/arion/projects/CommonMind/hoffman/sceqtl/sumStats/
R

library(ggplot2)
library(tidyverse)

assay_order = readRDS("assay_order.RDS")

ord.class = assay_order %>%
		grep("^class", ., value=TRUE) %>%
		gsub("^class / ", '', .)

ord.subclass = assay_order %>%
		grep("^subclass", ., value=TRUE) %>%
		gsub("^subclass / ", '', .)


# files = system("ls /sc/arion/projects/CommonMind/zengb02/single_cell_eQTL_for_Gabriel_02212024/eQTL_detection_on_pearsonsubtype_level/cis-eQTL_detection/analysis/colocalization_analysis/analysis/merged_colocalization_*_v2", intern=TRUE)

# df_coloc = lapply(files, function(file){

# 	ID = gsub("^.*merged_colocalization_(.*)$", "\\1", file )

# 	df = read_tsv(file)
# 	df$Trait = ID
# 	df[is.na(df)] = 0
# 	df
# 	})
# df_coloc = bind_rows(df_coloc) %>%
# 			pivot_longer(cols=!Gene&!Trait, names_to = "class", values_to="prob")

# fig = df_coloc %>%  
# 	filter(Trait == "alzBellenguez_v2",) %>%
# 	ggplot(aes(Gene, class, color = prob)) +
# 		geom_point() +
# 		theme_classic() +
# 		coord_equal() +
# 		scale_color_gradient(low="grey98", high="red1", limits=c(0,1)) +
# 		scale_x_discrete(guide = guide_axis(angle = 90)) +
# 		ggtitle("COLOC") +
# 		facet_wrap(~ Trait, ncol=1)

# ggsave(fig, file="plots/coloc_all.pdf", height = 10, width=9)




# Complete results
##################

files = system("ls /sc/arion/projects/CommonMind/zengb02/single_cell_eQTL_for_Gabriel_02212024/eQTL_detection_on_pearson*_level/cis-eQTL_detection/analysis/colocalization_analysis/results/merged_colocalization", inter=TRUE)

df = read.table(files[1], header=FALSE) %>%
		as_tibble 
colnames(df) = c("Trait", "CellType", "Gene", 'ppH0', 'ppH1', 'ppH2', 'ppH3', 'ppH4')

df = df %>%
		filter(Trait != "AD") 

write_tsv(df, file="coloc_class.tsv")


df_genes = df %>%
			group_by(Trait, Gene) %>%
			summarize(maxSignal = max(ppH4, na.rm=TRUE)) %>%
			filter(maxSignal > .5)

df_sub = df %>% 
		inner_join(df_genes, by=c("Trait", "Gene")) %>%
		select(Trait, CellType, Gene, ppH4)

res = lapply( unique(df_sub$Trait), function(x){
	M = df_sub %>%
		filter(Trait == x) %>%
		pivot_wider(names_from = "Gene", values_from=ppH4) %>%
		select(-Trait) %>% 
		column_to_rownames('CellType')
	M[is.na(M)] = 0

	hcl = hclust(dist(t(M)))

	gene.ord = hcl$labels[hcl$order]
})
gene.ord = unique(unlist(res))


fig = df_sub %>%  	
	mutate(CellType = factor(CellType, ord.class)) %>%
	mutate(Gene = factor(Gene, gene.ord)) %>%
	ggplot(aes(Gene, CellType, color = ppH4, size=ppH4, label=ifelse(ppH4 >=0.8, 'x', ''))) +
		geom_point() +
		theme_classic() +
		scale_color_gradient(low="grey98", high="red1", limits=c(0,1)) +
		scale_x_discrete(guide = guide_axis(angle = 90)) +
		ggtitle("COLOC") +
		geom_text(color="black", vjust=0.5, hjust=0.5) +
		scale_size_area(max_size=4) + 
		facet_wrap(~ Trait, nrow=1, scales="free_x")

file = paste0("plots/coloc_class.pdf")
ggsave(fig, file=file, height = 3.2, width=80, limitsize=FALSE)


# Subclass
##########

df = read.table(files[2], header=FALSE) %>%
		as_tibble 
colnames(df) = c("Trait", "CellType", "Gene", 'ppH0', 'ppH1', 'ppH2', 'ppH3', 'ppH4')

df = df %>%
		filter(Trait != "AD") 

write_tsv(df, file="coloc_subclass.tsv")

df_genes = df %>%
			filter(Trait != "AD") %>%
			group_by(Trait, Gene) %>%
			summarize(maxSignal = max(ppH4, na.rm=TRUE)) %>%
			filter(maxSignal > .5)

df_sub = df %>% 
		inner_join(df_genes, by=c("Trait", "Gene"))%>%
		select(Trait, CellType, Gene, ppH4)

res = lapply( unique(df_sub$Trait), function(x){
	M = df_sub %>%
		filter(Trait == x) %>%
		pivot_wider(names_from = "Gene", values_from=ppH4) %>%
		select(-Trait) %>% 
		column_to_rownames('CellType')
	M[is.na(M)] = 0

	hcl = hclust(dist(t(M)))

	gene.ord = hcl$labels[hcl$order]
})
gene.ord = unique(unlist(res))

fig = df_sub %>%  	
	mutate(CellType = factor(CellType, ord.subclass)) %>%
	mutate(Gene = factor(Gene, gene.ord)) %>%
	ggplot(aes(Gene, CellType, color = ppH4, size=ppH4, label=ifelse(ppH4 >=0.8, 'x', ''))) +
		geom_point() +
		theme_classic() +
		scale_color_gradient(low="grey98", high="red1", limits=c(0,1)) +
		scale_x_discrete(guide = guide_axis(angle = 90)) +
		ggtitle("COLOC") +
		geom_text(color="black", vjust=0.5, hjust=0.5) +
		scale_size_area(max_size=4) + 
		facet_wrap(~ Trait, nrow=1, scales="free_x")

file = paste0("plots/coloc_subclass.pdf")
ggsave(fig, file=file, height = 6, width=105, limitsize=FALSE)











# for( trait in unique(df_genes$Trait)){

# 	df_sub = df %>% 
# 		inner_join(df_genes, by=c("Trait", "Gene"))

# 	fig = df_sub %>%  	
# 		mutate(CellType = factor(CellType, ord.subclass)) %>%
# 		filter(Trait == trait) %>%
# 		ggplot(aes(Gene, CellType, color = ppH4, size=ppH4)) +
# 			geom_point() +
# 			theme_classic() +
# 			coord_equal() +
# 			scale_color_gradient(low="grey98", high="red1", limits=c(0,1)) +
# 			scale_x_discrete(guide = guide_axis(angle = 90)) +
# 			ggtitle("COLOC") +
# 			scale_size_area(max_size=4) + 
# 			facet_wrap(~ Trait, ncol=1)

# 	file = paste0("plots/coloc_", trait,".pdf")

# 	ggsave(fig, file=file, height = 18, width=8)
# }



