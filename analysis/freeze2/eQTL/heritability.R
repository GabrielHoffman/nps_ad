

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


# files = system("ls /sc/arion/projects/CommonMind/zengb02/single_cell_eQTL_for_Gabriel_02212024/eQTL_detection_on_pearson*level/heritability_analysis/estimate_heritability_explained/merged_expression_heritability*_level", inter=TRUE)


# df = read_tsv(files[2])

# fig = df %>% 
# 	ggplot(aes(cell, cis_heri, fill=cell)) +
# 		geom_violin() +
# 		theme_classic() +
# 		theme(aspect.ratio=2, legend.position="none")

# ggsave(fig, file="~/www/test.png")

# LDSC
#######

files = system("ls /sc/arion/projects/CommonMind/zengb02/single_cell_eQTL_for_Gabriel_02212024/eQTL_detection_on_pearson*_level/cis-eQTL_detection/analysis/finemapping_analysis/calculate_ld_score/run_ldsr/merged_LDSC_results", intern=TRUE)

df = read.table(files[2], header=FALSE) %>%
		as_tibble 
colnames(df) = c("CellType", "Trait", "effect", "se")

fig = df %>%
		filter(Trait != "ALS", Trait != "ALZ") %>%
		mutate(z = effect / se)  %>%
		arrange(-effect) %>%
		# filter(Trait %in% traits) %>%
		mutate(p.value = pnorm(z, lower.tail=FALSE)) %>%
		mutate(FDR = p.adjust(p.value)) %>%
		mutate(CellType = factor(CellType, ord.subclass)) %>%
		mutate(z = pmax(z, 0)) %>%
		mutate(effect = pmax(effect, 0)) %>%
	ggplot(aes(Trait, CellType, color = z, alpha = z, size = log10(effect+1), label=ifelse(FDR < 0.05, "x", ''))) +
		geom_point() +
		theme_classic() +
		coord_equal() +
		scale_color_gradient2(low="white", mid="grey98", high="red1") +
		scale_x_discrete(guide = guide_axis(angle = 90)) +
		scale_size_area() +
		geom_text(color="black", vjust=0.4, hjust=0.5, alpha = 1) 
ggsave(fig, file="plots/LDSC_subclass.pdf")


df = read.table(files[1], header=FALSE) %>%
		as_tibble 
colnames(df) = c("CellType", "Trait", "effect", "se")

fig = df %>%
		filter(Trait != "ALS", Trait != "ALZ") %>%
		mutate(z = effect / se)  %>%
		arrange(-effect) %>%
		# filter(Trait %in% traits) %>%
		mutate(p.value = pnorm(z, lower.tail=FALSE)) %>%
		mutate(FDR = p.adjust(p.value)) %>%
		mutate(CellType = factor(CellType, ord.class)) %>%
		mutate(z = pmax(z, 0)) %>%
		mutate(effect = pmax(effect, 0)) %>%
	ggplot(aes(Trait, CellType, color = z, alpha = z, size = log10(effect+1), label=ifelse(FDR < 0.05, "x", ''))) +
		geom_point() +
		theme_classic() +
		coord_equal() +
		scale_color_gradient2(low="white", mid="grey98", high="red1") +
		scale_x_discrete(guide = guide_axis(angle = 90)) +
		scale_size_area() +
		geom_text(color="black", vjust=0.4, hjust=0.5, alpha = 1) 
ggsave(fig, file="plots/LDSC_class.pdf")




# MESC
#######

files = system("ls /sc/arion/projects/CommonMind/zengb02/single_cell_eQTL_for_Gabriel_02212024/eQTL_detection_on_pearson*_level/MESC_analysis/run_MESC/merged_MESC_results", intern=TRUE)

files[2] = "/sc/arion/projects/CommonMind/zengb02/single_cell_eQTL_for_Gabriel_02212024/eQTL_detection_on_pearsonsubclass_level/MESC_analysis/run_MESC/another/merged_MESC_results"

df = read.table(files[2], header=FALSE) %>%
		as_tibble 
colnames(df) = c('Trait', 'CellType','h2_medi','h2_mediated_se', 'h2_percentage_mediated','h2_percentage_mediated_se')

fig = df %>% 
		filter(Trait != "ALS", Trait != "ALZ") %>%
		mutate(h2_percentage_mediated = pmax(0, h2_percentage_mediated)) %>%
		mutate(z = h2_percentage_mediated / h2_percentage_mediated_se)  %>%
		# filter(Trait %in% traits) %>%
		mutate(p.value = pnorm(z, lower.tail=FALSE)) %>%
		mutate(FDR = p.adjust(p.value)) %>%
		mutate(CellType = factor(CellType, ord.subclass)) %>% 
	ggplot(aes(Trait, CellType, color = z, alpha = z, size = h2_percentage_mediated, label=ifelse(FDR < 0.05, "x", ''))) +
		geom_point() +
		theme_classic() +
		coord_equal() +
		scale_color_gradient(low="grey98", high="red1") +
		scale_x_discrete(guide = guide_axis(angle = 90)) +
		scale_size_area(name="Hsq % mediated") +
		geom_text(color="black", vjust=0.4, hjust=0.5) +
		ggtitle("MESC")

ggsave(fig, file="plots/MESC_subclass.pdf", height = 7, width=7)


df = read.table(files[1], header=FALSE) %>%
		as_tibble 
colnames(df) = c('Trait', 'CellType','h2_medi','h2_mediated_se', 'h2_percentage_mediated','h2_percentage_mediated_se')

fig = df %>% 
		filter(Trait != "ALS", Trait != "ALZ") %>%
		mutate(h2_percentage_mediated = pmax(0, h2_percentage_mediated)) %>%
		mutate(z = h2_percentage_mediated / h2_percentage_mediated_se)  %>%
		# filter(Trait %in% traits) %>%
		mutate(p.value = pnorm(z, lower.tail=FALSE)) %>%
		mutate(FDR = p.adjust(p.value)) %>%
		mutate(CellType = factor(CellType, ord.class)) %>% 
	ggplot(aes(Trait, CellType, color = z, alpha = z, size = h2_percentage_mediated, label=ifelse(FDR < 0.05, "x", ''))) +
		geom_point() +
		theme_classic() +
		coord_equal() +
		scale_color_gradient(low="grey98", high="red1") +
		scale_x_discrete(guide = guide_axis(angle = 90)) +
		scale_size_area(name="Hsq % mediated") +
		geom_text(color="black", vjust=0.4, hjust=0.5) +
		ggtitle("MESC")

ggsave(fig, file="plots/MESC_class.pdf", height = 7, width=7)






