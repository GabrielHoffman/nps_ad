
# July 10, 2024
# adaped from nps_ad/analysis/freeze2/preprocess/preprocess.Rmd
# 
# Plot info about reads and cell count for each cell type


suppressPackageStartupMessages({
library(SingleCellExperiment)
library(zellkonverter)
library(DelayedArray)
library(tidyverse)
library(ggplot2)
library(cowplot) 
})



# order of cell types used for plotting
# Compute pseudobulk at multiple levels
ctorder = c("EN", 'EN_L2_3_IT', 'EN_L3_5_IT_1', 'EN_L3_5_IT_2', 'EN_L3_5_IT_3', 'EN_L5_6_NP', 'EN_L5_ET', 'EN_L6B', 'EN_L6_CT', 'EN_L6_IT_1', 'EN_L6_IT_2', "IN", 'IN_ADARB2', 'IN_LAMP5_LHX6', 'IN_LAMP5_RELN', 'IN_PVALB', 'IN_PVALB_CHC', 'IN_SST', 'IN_VIP', 'Astro', 'OPC', 'Oligo', 'Immune', 'Micro', 'PVM', 'Adaptive', 'Mural' ,'VLMC', 'SMC', 'PC', 'Endo')





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
names(cols) = gsub(".* / ", "", names(cols))



# read H5AD
path = "/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/h5ad_final/"
file = c(FULL = paste0(path, "FULL_2024-02-01_18_49.h5ad"))
 
sce = readH5AD(file, use_hdf5=TRUE, verbose=TRUE, version="0.8.0")




# Number of cells observed per Subject
figA = colData(sce) %>%
  xtabs( ~ SubID + class,.) %>%
  as_tibble %>%
  pivot_longer(cols=SubID) %>%
  mutate(class = droplevels(factor(class, ctorder))) %>%  
  ggplot(aes(class, n+.1, fill=class)) + 
    geom_violin(color = NA) + 
    geom_boxplot(width=.1, outlier.size=.1) +
    theme_bw() + 
    theme(aspect.ratio=1.5, plot.title = element_text(hjust = 0.5),
      legend.position="none") + 
    scale_y_log10(limits=c(1, 10000), expand=c(0,0)) +
    coord_flip() +
    ylab("Number of cells observed per Subject") +
    xlab('') +
    scale_fill_manual(values = cols)

figB = colData(sce) %>%
  xtabs( ~ SubID + subclass,.) %>%
  as_tibble %>%
  pivot_longer(cols=SubID) %>%
  mutate(subclass = factor(subclass, ctorder)) %>%  
  ggplot(aes(subclass, n+.1, fill=subclass)) + 
    geom_violin(color = NA) + 
    geom_boxplot(width=.1, outlier.size=.1) +
    theme_bw() + 
    theme(aspect.ratio=1.5, plot.title = element_text(hjust = 0.5),
      legend.position="none") + 
    scale_y_log10(limits=c(1, 10000), expand=c(0,0)) +
    coord_flip() +
    ylab("Number of cells observed per Subject") +
    xlab('') +
    scale_fill_manual(values = cols)

# Number of cells observed per Subject
figC = colData(sce) %>%
  xtabs( ~ SubID + subtype,.) %>%
  as_tibble %>%
  pivot_longer(cols=SubID) %>%
  # mutate(subtype = factor(subtype, ctorder)) %>%  
  ggplot(aes(subtype, n+.1, fill=subtype)) + 
    geom_violin(color = NA) + 
    geom_boxplot(width=.1, outlier.size=.1) +
    theme_bw() + 
    theme(aspect.ratio=3, plot.title = element_text(hjust = 0.5),
      legend.position="none") + 
    scale_y_log10(limits=c(1, 10000), expand=c(0,0)) +
    coord_flip() +
    ylab("Number of cells observed per Subject") +
    xlab('') +
    scale_fill_manual(values = cols)

# fig = plot_grid(figA, figB, nrow=1, align="v")
# ggsave(fig, file="~/www/Cells_per_subject_AB.pdf", width=10)


ggsave(figC, file="~/www/Cells_per_subject_C.pdf", width=15, height=10)



library(dreamlet)

# Read from RDS

files = dir("/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/pseudobulk/", pattern="FULL_2024-02-.*_PB_SubID_.*.RDS", full.names=TRUE)

# reads per cell
res = lapply(files[-1], function(file){

	pb = readRDS(file)

	# extract read counts for each Channel
	df_counts = lapply( assayNames(pb), function(x){
	  data = assay(pb, x)

	  data.frame(celltype = x, ID = colnames(data), readCounts = colSums(data))
	})
	df_counts = do.call(rbind, df_counts)
	df_counts$celltype = factor(df_counts$celltype, ctorder)  

	# extract cell counts
	df_rate = cellCounts(pb) %>%
	            as.data.frame %>%
	            mutate(ID = rownames(.))  %>% 
	            pivot_longer(cols=-ID, values_to="ncells", names_to="celltype") %>%
	            mutate(celltype = factor(celltype, ctorder))  

	# plot reads per cell
	inner_join(df_counts, df_rate, by=c("celltype", "ID")) %>%
	  mutate(celltype = droplevels(factor(celltype, ctorder))) %>%
	  ggplot(aes(celltype, readCounts/ncells, fill=celltype)) +  
	    geom_violin(color = NA) + 
	    geom_boxplot(width=.1, outlier.size=.1) +
	    theme_bw() + 
	    theme(aspect.ratio=3, plot.title = element_text(hjust = 0.5),
	      legend.position="none") + 
	    scale_y_log10() +
	    coord_flip() +
	    ylab('Reads per cell') +
	    xlab('') +
	    ggtitle('Reads per cell') +
	    scale_fill_manual(values = cols)
})

fig = plot_grid(plotlist = res[1:2], nrow=1, align="v")
ggsave(fig, file="~/www/read_per_cell_AB.pdf", width=15, height=10)

ggsave(res[[3]], file="~/www/read_per_cell_C.pdf", width=15, height=10)



# reads per Subject
res = lapply(files[-1], function(file){

	pb = readRDS(file)

	# extract read counts for each Channel
	df_counts = lapply( assayNames(pb), function(x){
	  data = assay(pb, x)

	  data.frame(celltype = x, ID = colnames(data), readCounts = colSums(data))
	})
	df_counts = do.call(rbind, df_counts)
	df_counts$celltype = factor(df_counts$celltype, ctorder)  

	# extract cell counts
	df_rate = cellCounts(pb) %>%
	            as.data.frame %>%
	            mutate(ID = rownames(.))  %>% 
	            pivot_longer(cols=-ID, values_to="ncells", names_to="celltype") %>%
	            mutate(celltype = factor(celltype, ctorder))  

	# plot reads per cell
	inner_join(df_counts, df_rate, by=c("celltype", "ID")) %>%
	  mutate(celltype = droplevels(factor(celltype, ctorder))) %>%
	  ggplot(aes(celltype, readCounts, fill=celltype)) +  
	    geom_violin(color = NA) + 
	    geom_boxplot(width=.1, outlier.size=.1) +
	    theme_bw() + 
	    theme(aspect.ratio=3, plot.title = element_text(hjust = 0.5),
	      legend.position="none") + 
	    scale_y_log10() +
	    coord_flip() +
	    ylab('Reads per cell') +
	    xlab('') +
	    ggtitle('Reads per cell') +
	    scale_fill_manual(values = cols)
})

fig = plot_grid(plotlist = res[1:2], nrow=1, align="v")
ggsave(fig, file="~/www/read_per_Subject_AB.pdf", width=15, height=10)

ggsave(res[[3]], file="~/www/read_per_Subject_C.pdf", width=15, height=10)


# Save mean reads per Subject
df_counts = lapply(files, function(file){

	pb = readRDS(file)

	# extract read counts for each Channel
	df_counts = lapply( assayNames(pb), function(x){
	  data = assay(pb, x)

	  tibble(celltype = x, ID = colnames(data), readCounts = colSums(data))
	})
	df_counts = bind_rows(df_counts) %>%
			mutate(Level = gsub("^.*_(\\S+).RDS", "\\1", basename(file)))
	
	df_counts %>%
		group_by(celltype, Level) %>%
		summarize(meanReadCount = mean(readCounts))	 
})
df_counts = bind_rows(df_counts)

outFile = "/sc/arion/projects/CommonMind/hoffman/sceqtl/sumStats/df_counts.RDS"
saveRDS(df_counts, file=outFile)


