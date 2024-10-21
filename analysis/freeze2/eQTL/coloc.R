



cd /sc/arion/projects/CommonMind/hoffman/sceqtl/sumStats/
R

library(ggplot2)
library(tidyverse)
library(ggh4x)

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


table(df$Trait)

df = df %>%
		filter(Trait != "AD") %>%
		mutate(Trait = factor(Trait)) %>%
		mutate(Trait = fct_recode(Trait, AD = 'AD3', SCZ = 'SCZ3', MDD = "MDD2"))

write_tsv(df, file="coloc_class.tsv")



height = list('0.5' = 40, '0.8' = 13)
width = list('0.5' = 3.5, '0.8' = 3.2)

res = lapply( c(0.5, 0.8), function(cutoff){
	df_genes = df %>%
				group_by(Trait, Gene) %>%
				summarize(maxSignal = max(ppH4, na.rm=TRUE)) %>%
				filter(maxSignal > cutoff)

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

	df_sub2 = df_sub %>%  	
		mutate(CellType = factor(CellType, ord.class)) %>%
		mutate(Gene = factor(Gene, gene.ord)) %>%	
		mutate(label=ifelse(ppH4 >=0.5, 'o', '')) %>%
		mutate(label=ifelse(ppH4 >=0.8, 'x', label)) 

	counts = df_sub2 %>%
		group_by(Trait) %>%
		select(Gene) %>%
		distinct %>%
		summarize(n = n())

	fig = df_sub %>%  	
		mutate(CellType = factor(CellType, ord.class)) %>%
		mutate(Gene = factor(Gene, gene.ord)) %>%
		mutate(label=ifelse(ppH4 >=0.5, '*', '')) %>%
		mutate(label=ifelse(ppH4 >=0.8, 'x', label)) %>%
		ggplot(aes(CellType, Gene, color = ppH4, size=ppH4, label=label)) +
			geom_point() +
			theme_classic() +
			scale_color_gradient(low="grey98", high="red1", limits=c(0,1)) +
			scale_x_discrete(guide = guide_axis(angle = 90)) +
			ggtitle("COLOC") +
			geom_text(color="black", vjust=0.5, hjust=0.5) +
			scale_size_area(max_size=4) +  
			facet_wrap(~ Trait, ncol=1, scales="free_y") +
			force_panelsizes(rows = counts$n)  +
			theme(plot.title = element_text(hjust = 0.5))

	file = paste0("plots/coloc_class_", cutoff, ".pdf")
	ggsave(fig, file=file, height = height[[as.character(cutoff)]], width=width[[as.character(cutoff)]], limitsize=FALSE)

	library(cowplot)

	fig1 = df_sub %>%
		group_by(Trait, CellType) %>%
		select(Gene) %>%
		distinct  %>%
		summarize(n = n()) %>%
		ggplot(aes(n, CellType, fill = CellType)) +
			geom_bar(stat="identity") +
			facet_wrap(~Trait, nrow=1) +
			theme_classic() +
			theme(aspect.ratio=1, legend.position="none") +
			scale_x_continuous(limits=c(0, NA), expand=c(0,0))

	fig2 = df_sub %>%
		group_by(Trait) %>%
		select(Gene) %>%
		distinct  %>%
		summarize(n = n()) %>%
		ggplot(aes(n, Trait, fill = Trait)) +
			geom_bar(stat="identity") +
			theme_classic() +
			theme(aspect.ratio=1, legend.position="none")+
			scale_x_continuous(limits=c(0, NA), expand=c(0,0))

	fig = plot_grid(fig1, fig2, rel_widths = c(1,.3))
	file = paste0("plots/coloc_class_count_", cutoff, ".pdf")
	ggsave(fig, file=file, height = 4, width=7, limitsize=FALSE)
})









# Subclass
##########

df = read.table(files[2], header=FALSE) %>%
		as_tibble 
colnames(df) = c("Trait", "CellType", "Gene", 'ppH0', 'ppH1', 'ppH2', 'ppH3', 'ppH4')

df = df %>%
		filter(Trait != "AD") %>%
		mutate(Trait = factor(Trait)) %>%
		mutate(Trait = fct_recode(Trait, AD = 'AD3', SCZ = 'SCZ3', MDD = "MDD2"))

write_tsv(df, file="coloc_subclass.tsv")

height = list('0.5' = 52, '0.8' = 18)

res = lapply( c(0.5, 0.8), function(cutoff){
	df_genes = df %>%
				group_by(Trait, Gene) %>%
				summarize(maxSignal = max(ppH4, na.rm=TRUE)) %>%
				filter(maxSignal > cutoff)

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


	# All colocs
	#-----------
	df_sub2 = df_sub %>%  	
		mutate(CellType = factor(CellType, ord.subclass)) %>%
		mutate(Gene = factor(Gene, gene.ord)) %>%	
		mutate(label=ifelse(ppH4 >=0.5, 'o', '')) %>%
		mutate(label=ifelse(ppH4 >=0.8, 'x', label)) 

	counts = df_sub2 %>%
		group_by(Trait) %>%
		select(Gene) %>%
		distinct %>%
		summarize(n = n())

	fig = df_sub2 %>%
		ggplot(aes(CellType, Gene, color = ppH4, size=ppH4, label=label)) +
			geom_point() +
			theme_classic() +
			scale_color_gradient(low="grey98", high="red1", limits=c(0,1)) +
			scale_x_discrete(guide = guide_axis(angle = 90)) +
			ggtitle("COLOC") +
			geom_text(color="black", vjust=0.4, hjust=0.4) +
			scale_size_area(max_size=4) +  
			facet_wrap(~ Trait, ncol=1, scales="free_y") +
			force_panelsizes(rows = counts$n) +
			theme(plot.title = element_text(hjust = 0.5))

	file = paste0("plots/coloc_subclass_", cutoff, ".pdf")
	ggsave(fig, file=file, height = height[[as.character(cutoff)]], width=5.8, limitsize=FALSE)

	fig1 = df_sub %>%
		group_by(Trait, CellType) %>%
		select(Gene) %>%
		distinct  %>%
		summarize(n = n()) %>%
		ggplot(aes(n, CellType, fill = CellType)) +
			geom_bar(stat="identity") +
			facet_wrap(~Trait, nrow=1) +
			theme_classic() +
			theme(aspect.ratio=1, legend.position="none") +
			scale_x_continuous(limits=c(0, NA), expand=c(0,0))

	fig2 = df_sub %>%
		group_by(Trait) %>%
		select(Gene) %>%
		distinct  %>%
		summarize(n = n()) %>%
		ggplot(aes(n, Trait, fill = Trait)) +
			geom_bar(stat="identity") +
			theme_classic() +
			theme(aspect.ratio=1, legend.position="none")+
			scale_x_continuous(limits=c(0, NA), expand=c(0,0))

	fig = plot_grid(fig1, fig2, rel_widths = c(1,.2))
	file = paste0("plots/coloc_subclass_count_", cutoff, ".pdf")
	ggsave(fig, file=file, height = 4, width=15, limitsize=FALSE)
})



# Compare Class vs Subclass
###########################

# Genes that coloc in Subclass but not class
df_class = read_tsv("coloc_class.tsv", show_col_types = FALSE) %>%
			group_by(Trait, Gene) %>%
			summarize(maxSignal = max(ppH4, na.rm=TRUE)) %>%
			filter(maxSignal > 0.8)

df_subclass = read_tsv("coloc_subclass.tsv", show_col_types = FALSE)

res = lapply( unique(df_subclass$Trait), function(TR){

	 df_subclass %>%
		filter(Trait == TR, ppH4 > 0.8) %>%
		select(-ppH0, -ppH1, -ppH2, -ppH3, -ppH4, -CellType) %>%left_join( df_class) %>%
		distinct %>%
		arrange(maxSignal) %>%
		filter(is.na(maxSignal))
})
res = bind_rows(res)


res %>% 
	write_tsv(file="coloc_subclass_not_class.tsv")


res %>%
	group_by(Trait) %>%
	summarize(n = n())


# Compare class and subclass overlap
df1 = read_tsv("coloc_class.tsv", show_col_types = FALSE) %>%
	filter(ppH4 > 0.8) %>%
	group_by(Trait) %>%
	select(Gene) %>%
	distinct %>%
	mutate(Discovered = "class") 

df2 = df_subclass %>%
	filter(ppH4 > 0.8) %>%
	group_by(Trait) %>%
	select(Gene) %>%
	distinct %>%
	mutate(Discovered = "subclass") 

df = bind_rows(df1, df2)


res = lapply( unique(df$Trait), function(TR){
	df1 = df %>%
		filter(Trait == TR, Discovered == "class") 

	df2 = df %>%
		filter(Trait == TR, Discovered == "subclass") 

	tibble(Trait = TR,
			shared = length(intersect(df1$Gene, df2$Gene)),
			class = length(setdiff(df1$Gene, df2$Gene)),
			subclass = length(setdiff(df2$Gene, df1$Gene)))
})
res = bind_rows(res)


cols = RColorBrewer::brewer.pal(8, "Set1")
names(cols)[1] = "AD"
names(cols)[2] = "ASD"
names(cols)[3] = "MDD"
names(cols)[4] = "SCZ"


fig = res %>%
	pivot_longer(!Trait) %>%
	mutate(name = fct_relevel(name, c("subclass", "class", "shared"))) %>%
	ggplot(aes(Trait, value, fill=Trait, alpha=name)) +
		geom_bar(stat="identity") +
		theme_classic() +
		theme(aspect.ratio=1) +
		ylab("# colocalized genes") +
		scale_fill_manual(values=cols[1:4]) +
		scale_y_continuous(expand=c(0,0))

file = "plots/coloc_overlap.pdf"
ggsave(fig, file=file, height = 6, width=6) 












	%>%
	with(table(Gene, Discovered)) %>%
	data.frame


fig = df %>%

		summarize(n )


ggplot(df, aes(Trait))






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



