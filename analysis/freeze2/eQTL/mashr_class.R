

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


file = "/sc/arion/projects/CommonMind/zengb02/single_cell_eQTL_for_Gabriel_02212024/eQTL_detection_on_pearsonclass_level/cis-eQTL_detection/analysis/evaluate_eQTL_sharing_mashR/get_lfsr/another2/lfsr_results_for_psychAD_sceQTL_class"

# CLASS
df_class = read.table(file, row.names=1, header=TRUE) %>%
				rownames_to_column("Gene") %>%
				as_tibble
colnames(df_class) = colnames(df_class) %>% 
						gsub("_neurons", '', .) %>% 
						gsub("EX", "EN", .) %>% 
						gsub("Microglia", "Immune", .) %>% 
						gsub("Astrocytes", "Astro", .) %>% 
						gsub("Oligodendrocyte", "Oligo", .) 

df_class %>%
	write_tsv(file="mashr_lfsr_class.tsv.gz")

res = apply(df_class[,-1], 2, function(x) which(x < 0.05))

ord = assay_order %>%
		grep("^class", ., value=TRUE) %>%
		gsub("^class / ", '', .)

# ord = c("IN", "EN", "Oligo", "OPC", "Astro", "Immune")
ord = assay_order %>%
		grep("^class", .,value=TRUE) %>%
		gsub("class / ", "", .)

pdf("plots/upset_mashr_class.pdf", height=5, width=7)
upset(fromList(res[ord]), order.by = "freq", nsets=8, keep.order=TRUE, sets=ord, nintersects=15)
dev.off()






# Strong cell type specificity
##############################
file = c("/sc/arion/projects/CommonMind/zengb02/single_cell_eQTL_for_Gabriel_02212024/eQTL_detection_on_pearsonclass_level/cis-eQTL_detection/analysis/evaluate_eQTL_sharing_mashR/get_lfsr/another2/lfsr_results_for_psychAD_sceQTL_class")

# read in results
df_lfsr = read.table(file, row.names=1, header=TRUE) %>%
				rownames_to_column("Gene") %>%
				as_tibble %>%
		mutate(ID = Gene) %>%
		mutate(Gene = gsub("^(\\S+)_.*$", "\\1", Gene))



library(dreamlet)
ctlist = as.list(colnames(df_lfsr)[c(-1,-length(df_lfsr))])
ctlist[[length(ctlist)+1]] = c("EN", "IN")

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

# number of unique genes
a = "specific_eQTL_class.tsv.gz" %>%
	read_tsv %>%
	filter(prob > 0.5) %>%
	select(Gene) %>%
	distinct 
nrow(a)

# number of unique genes
b = "specific_eQTL_subclass.tsv.gz" %>%
	read_tsv %>%
	filter(prob > 0.5) %>%
	select(Gene) %>%
	distinct 

nrow(b)

c(a$Gene, gsub("^(\\S+)_.*$", "\\1", b$Gene)) %>%
	unique %>%
	length

"specific_eQTL_subclass.tsv.gz" %>%
	read_tsv %>%
	filter(Gene == "EGFR_rs74504435")



"specific_eQTL_subclass.tsv.gz" %>%
	read_tsv %>%
	filter(grepl("EGFR_", Gene)) %>%
	arrange(-prob)








