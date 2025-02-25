


cd /sc/arion/projects/CommonMind/hoffman/sceqtl/sumStats/
R


library(ggplot2)
library(tidyverse)
library(UpSetR)
library(arrow)

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


#################################
# Get ALL trans-eQTL to parquet #
#################################

library(ggplot2)
library(tidyverse)
library(arrow)
library(parallel)

# Write trans eQTL results to parquet
files = system("ls /sc/arion/projects/CommonMind/zengb02/single_cell_eQTL_for_Gabriel_02212024/eQTL_detection_on_pearsonclass_level/trans-eQTL_detection/another/eQTL_results/eQTL_result_MSSM_*", intern=TRUE)

df = mclapply(files, function(file){
	
	dfin = readr::read_tsv(file, col_names = FALSE, show_col_types=FALSE, progress=FALSE, num_threads=4)
	colnames(dfin) = c("chrom", 'start', 'end', 'variant', 'allele1', 'allele2_assessed', 'gene', 'beta', 'se', 'z')

	CellType = gsub("eQTL_result_MSSM_(\\S+)_.*_.*$", "\\1", basename(file))
	dfin$CellType = CellType
	
	chrom = paste0("chr", dfin$chrom[1])

	i = match(file, files)
	outfile = paste0("tmp/trans_eQTL_results_", chrom, "_", sprintf("%04d", i),".parquet")
	arrow::write_parquet(dfin, outfile, compression="UNCOMPRESSED")
	1
	}, mc.cores=12)


files = system("ls tmp/trans_eQTL_results_*.parquet", intern=TRUE)

txt = paste0("_chr", seq(22), "_")
res = mclapply(txt, function(chrom){

	i = grep(chrom, files)

	dSet = open_dataset(files[i])

	outfile = paste0("export/trans/trans_eQTL_results", chrom, ".parquet")
	write_parquet(dSet, outfile)
	colnames(dSet)
}, mc.cores=12)




############################
# Select Significant Genes #
############################

# Get ENS -> HGNC mapping
library(biomaRt)
library(tidyverse)
mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
geneInfo = getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position"), mart = mart)   
write_tsv(geneInfo, file="geneInfo.tsv")  

cd /sc/arion/projects/CommonMind/hoffman/sceqtl/sumStats/
R

library(tidyverse)
library(arrow)

# Mappability
# file = "extern/hg38_cross_mappability_strength_symmetric_mean_sorted.txt.gz"
# df_map = read_tsv(file, col_names=c("Ens1", "Ens2", "Score"), show_col_types=FALSE) %>%
# 	mutate(Ens1 = gsub("^(\\S+)\\.(\\S+)$", "\\1", Ens1)) %>%
# 	mutate(Ens2 = gsub("^(\\S+)\\.(\\S+)$", "\\1", Ens2))


file = "extern/hg38_gene_mappability.txt.gz"
df_map = read_tsv(file, col_names=c("EnsID", "Score"), show_col_types=FALSE) %>%
	filter(!grepl("PAR_Y",EnsID)) %>%
	mutate(EnsID = gsub("^(\\S+)\\.(\\S+)$", "\\1", EnsID)) 

# file = "extern/hg38_snp_mappability_75mer.bed.gz"
# df_snp = read_tsv(file, col_names = c("chrom", "start", "end", "Score"), skip=1, n_max=100, show_col_types = FALSE) 


geneInfo = read_tsv("geneInfo.tsv", show_col_types = FALSE) %>%
				rename(EnsID = ensembl_gene_id, 
					Gene = hgnc_symbol, 
					chrom.gene = chromosome_name, 
					start.gene = start_position, 
					end.gene = end_position) %>%
				filter(EnsID != "ENSG00000291104") %>%
				filter(chrom.gene %in% c(1:22, "X", "Y", "MT")) %>%
				left_join(df_map)

# indicate if gene has map info
# ids = geneInfo$EnsID[geneInfo$EnsID %in% unique(df_map$Ens1)]
# geneInfo$mapInfo = FALSE
# geneInfo$mapInfo[geneInfo$EnsID %in% ids] = TRUE

write_tsv(geneInfo, file="geneInfo.tsv.gz")



# df = read_parquet("trans_eQTL_results.parquet")

files = system("ls export/trans/trans_eQTL_results*.parquet", intern=TRUE)
df = open_dataset(files)

#  2,716,381,199 tests

# format(nrow(df), big.mark=',')



# # cutoff 
# pnorm(5.0, lower.tail=F)

# keep SNP -> Gene with z > 5
# Bonferroni test
# qnorm(0.05/nrow(df), lower.tail=FALSE)
df2 = df %>%
		filter(abs(z) > 5) %>%
		collect

# Number of t
df_test = df %>% 
		rename(Gene = gene) %>% 
		group_by(CellType, Gene) %>%
		summarize(ntests = n()) %>%
		collect 


# a = df %>%
# 	filter(gene == 'SAMSN1' , CellType == 'Immune') %>%
# 	collect

# Get CellType, Variant, Gene
# for most significant Variant for CellType/Gene
ntests = 1e5
df_gene = df2 %>%
			rename(Variant = variant, Gene = gene) %>%
			group_by(CellType, Gene) %>%
			slice(which.max(abs(z))) %>%
			left_join(geneInfo %>%
				group_by(Gene) %>%
				slice(which.min(Score)), 
				by="Gene") %>%
			mutate(p.value = 2*pnorm(abs(z), lower.tail=FALSE)) %>%
			mutate(p.sidak = 1 - (1 - min(p.value, na.rm=TRUE))^ntests) %>%
			filter(Score > 0.8) %>%
			mutate(FDR = p.adjust(p.sidak, "BH")) %>%
			filter(FDR < 0.05) 


write_tsv(df_gene, file="export/trans/trans_eQTL_FDR05.tsv.gz")



# count trans-eGenes
#######################

cd /sc/arion/projects/CommonMind/hoffman/sceqtl/sumStats/
R

library(tidyverse)
library(ggplot2)

df_genes = read_tsv("export/trans/trans_eQTL_FDR05.tsv.gz")

fig = df_genes %>%
		group_by(CellType) %>%
		summarize(n = n()) %>%
		mutate(CellType = factor(CellType, ord)) %>%
		ggplot(aes(n, CellType, fill = CellType, label=n)) +
			geom_bar(stat="identity") +
			scale_x_continuous(limits=c(0, 420), expand=c(0,0)) +
			theme_classic() +
			theme(aspect.ratio=1, legend.position="none") +
			scale_fill_manual(values = cols) +
			geom_text() +
			xlab("# genes with trans-eQTL")

ggsave(fig, file="plots/trans_eQTL_count.pdf")


res = lapply( unique(df_genes$CellType), function(CT){
	df_genes %>%
		filter(CellType == CT) %>%
		pull(Gene) %>%
		unique
	})
names(res) = unique(df_genes$CellType)


ord = ord[ord %in% names(res)]
pdf("plots/trans_upset.pdf", height=5, width=5)
upset(fromList(res[ord]), order.by = "freq", nsets=8, keep.order=TRUE, sets=ord, nintersects=15)
dev.off()

# CPSF4L

lapply(1:22, function(chr){
	df_genes %>%
		filter(chrom == chr) %>%
		arrange(start) %>% 
		head(10) %>% data.frame
	})


get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


library(viridis)

fig = df_genes %>%
		mutate(density = get_density(chrom, start)) %>%
		ggplot(aes(start, abs(z), color=density)) +
			geom_point() +
			theme_classic() +
			theme(legend.position="none") +
			facet_grid(chrom ~ CellType) +
			ylim(4, NA) + scale_color_viridis()


ggsave(fig, file="plots/trans_location.pdf", height=20, width=12)



fig = geneInfo %>%
		ggplot(aes(start.gene, Score)) +
			geom_point() +
			theme_classic() +
			theme(legend.position="none") +
			facet_wrap(~ chrom.gene, ncol=3) +
			ylim(0, 1) + scale_color_viridis()


ggsave(fig, file="plots/geneInfo_location.pdf", height=20, width=12)


# variants with multiple targets
df_trans_multi = df_genes %>%
	group_by(CellType, Variant) %>%
	summarize(n = n()) %>%
	filter(n > 2) %>%
	arrange(-n)

# df_genes %>%
# 	filter(CellType == "Oligo", chrom == 6) %>%
# 	arrange(start) %>%
# 	select(start, Variant, Gene)  %>%
# 	data.frame


# CIRCOS plot
#############

library(RCircos)
data(UCSC.HG38.Human.CytoBandIdeogram)
GRCm38.cyto <- UCSC.HG38.Human.CytoBandIdeogram

for(i in seq(nrow(df_trans_multi))){
	CT = df_trans_multi$CellType[i]
	v = df_trans_multi$Variant[i]

	df_sub = df_genes %>%
				filter(CellType == CT, Variant == v) %>%
				mutate(chrom = paste0("chr", chrom))%>%
				mutate(chrom.gene = paste0("chr", chrom.gene)) %>%
		rename(	Chromosome = chrom, 
				chromStart = start, 
				chromEnd = end,
				Chromosome.1 = chrom.gene, 
				chromStart.1 = start.gene, 
				chromEnd.1 = end.gene) %>%
		select(Chromosome, chromStart, chromEnd,
			Chromosome.1, chromStart.1, chromEnd.1, Gene) %>%
		data.frame

	file = paste0("plots/trans_circos_", CT, "_", v, ".pdf")
	pdf(file)
	RCircos.Set.Core.Components(cyto.info=GRCm38.cyto, chr.exclude=c("chrX", "chrY"), tracks.inside=3, tracks.outside=0)
	RCircos.Set.Plot.Area()
	RCircos.Chromosome.Ideogram.Plot();
	RCircos.Link.Plot(link.data=df_sub, track.num=2, by.chromosome=FALSE);

	RCircos.Gene.Name.Plot(gene.data = df_sub[,4:7], 4, 2, "in")
	RCircos.Gene.Connector.Plot(df_sub[,4:7],1,"in")
	dev.off()
}










df3 = df_genes %>%
	filter(chrom == 17) %>%
	filter(start > 40160522, end < 49944823) 







df_mapt = df2 %>%
	rename(Variant = variant, Gene = gene) %>%
	filter(chrom == 17) %>%
	filter(start > 40160522, end < 49944823) %>%
	select(CellType, chrom, start, end, Variant, Gene, z ) %>%
	group_by(CellType) %>%
	reframe(Gene = unique(Gene))



df_mapt$Gene %in% unique(dfnew$gene)
unique(dfnew$gene)[!unique(dfnew$gene)%in%df_mapt$Gene]


"TXNL4B" %in% df_mapt$Gene


# All variants in MAPT
# /sc/arion/projects/CommonMind/zengb02/single_cell_eQTL_for_Gabriel_02212024/eQTL_detection_on_pearsonclass_level/trans-eQTL_detection/another/for_MAPT_region/eQTL_results/*.gz

# looked at all the variants in MAPT
files = system("ls /sc/arion/projects/CommonMind/zengb02/single_cell_eQTL_for_Gabriel_02212024/eQTL_detection_on_pearsonclass_level/trans-eQTL_detection/another/for_MAPT_region/analysis/significant_trans_eQTL_in_*", intern=TRUE)

cn = c("chrom", 'start', 'end', 'variant', 'allele1', 'allele2_assessed', 'Gene', 'beta', 'se', 'z')
df_mapt = lapply(files, function(file){
	CT = gsub("^.*_(\\S+)$", "\\1", basename(file))
	read_tsv(file, col_names = cn, show_col_types=FALSE) %>%
		mutate(CellType = CT)
	}) %>%
	bind_rows 


df_mapt %>%
	filter(variant == "rs17769552") %>%
		left_join(geneInfo %>%
				group_by(Gene) %>%
				slice(which.min(Score)), by="Gene") %>%
	filter(Score > 0.8) %>%
	data.frame

df_mapt = df_mapt %>%
			left_join(geneInfo %>%
				group_by(Gene) %>%
				slice(which.min(Score)), 
				by="Gene") %>%
			filter(Score > .8, chrom.gene != "chr17")



res = lapply( unique(df_mapt$CellType), function(CT){
	df_mapt %>%
		filter(CellType == CT) %>%
		pull(Gene) %>%
		unique
	})
names(res) = unique(df_mapt$CellType)


ord = ord[ord %in% names(res)]
pdf("plots/trans_upset_MAPT.pdf", height=5, width=5)
upset(fromList(res[ord]), order.by = "freq", nsets=8, keep.order=TRUE, sets=ord, nintersects=15)
dev.off()



# CIRCOS plot
#############


# variants with multiple targets
df_trans_multi = df_mapt %>%
	group_by(CellType, variant) %>%
	summarize(n = n()) %>%
	arrange(-n)


library(RCircos)
data(UCSC.HG38.Human.CytoBandIdeogram)
GRCm38.cyto <- UCSC.HG38.Human.CytoBandIdeogram

for(CT in unique(df_mapt$CellType)){

	df_sub = df_mapt %>%
				filter(CellType == CT) %>%
				mutate(chrom = paste0("chr", chrom))%>%
				mutate(chrom.gene = paste0("chr", chrom.gene)) %>%
		rename(	Chromosome = chrom, 
				chromStart = start, 
				chromEnd = end,
				Chromosome.1 = chrom.gene, 
				chromStart.1 = start.gene, 
				chromEnd.1 = end.gene) %>%
		select(Chromosome, chromStart, chromEnd,
			Chromosome.1, chromStart.1, chromEnd.1, Gene) %>%
		group_by(Gene) %>%
		summarize(Chromosome = Chromosome[1], chromStart = min(chromStart), chromEnd = max(chromEnd), Chromosome.1 = Chromosome.1[1],  chromStart.1=chromStart.1[1], chromEnd.1 = chromEnd.1[1]) %>%
		filter(Chromosome.1 != "chr17") %>%
		data.frame

	file = paste0("plots/trans_circos_MAPT_", CT,".pdf")
	pdf(file)
	RCircos.Set.Core.Components(cyto.info=GRCm38.cyto, chr.exclude=c("chrX", "chrY"), tracks.inside=3, tracks.outside=0)
	RCircos.Set.Plot.Area()
	RCircos.Chromosome.Ideogram.Plot();
	RCircos.Link.Plot(link.data=df_sub[,-1], track.num=2, by.chromosome=FALSE);

	RCircos.Gene.Name.Plot(gene.data = df_sub[,c(5:7,1)], 4, 2, "in")
	RCircos.Gene.Connector.Plot(df_sub[,c(5:7,1)],1,"in")
	dev.off()
}

























sort(unique(dfnew$gene))

a = df %>%
	filter(variant == "rs16966017") %>%
	collect




df_genes %>%
	select(Gene) %>%
	distinct %>%
			left_join(geneInfo, by="Gene")

df_genes %>%
	filter(Gene == "MAPT")

# MAPT region
#  17: 45,894,527-46,028,334 




# count trans-eGenes
df_trans_genes = df_genes %>%
	group_by(CellType) %>%
	reframe(Gene = unique(Gene), Variant) 


df_trans_genes %>%
	group_by(CellType) %>%
	summarize(n = length(unique(Gene)))










df_trans_genes = df_genes %>%
	select(Variant, Gene, z, CellType, EnsID, start.gene, end.gene, mapInfo)

a = df_trans_genes %>%
	group_by(CellType) %>%
	summarize(n = length(unique(Gene)))



a = df_trans_genes %>%
	group_by(CellType, Gene) %>%
	summarize(n = length(Gene)) %>%
	arrange(-n)



gene = "ADAM3A"

df_trans_genes %>%
	filter(Gene == gene)


df_genes %>%
	filter(Gene == gene) %>%
	data.frame




###########
# FIGURES #
###########

# Dropped MAPT from chr17, but can we type that CNV

cd /sc/arion/projects/CommonMind/hoffman/sceqtl/sumStats/
# cd /Users/gabrielhoffman/Library/CloudStorage/Dropbox/projects/PsychAD3.0/supplementary_figures/trans-eQTL
R

library(ggplot2)
library(tidyverse)
library(UpSetR)
library(arrow)

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

# df_all = open_dataset("trans_eQTL_results.parquet")

# Plot genes, or gene/variant pairs ?

df_top = read_tsv("trans_eQTL_top_hits.tsv") %>%
			rename(Gene = `trans-eGENE`, 
				Chrom = "chr_eGene", 
				Start = "eGene_start", 
				End = "eGene_start",
				Variant = `trans-eQTL`,
				ChromVariant = "chr_eQTL",
				StartVariant = "eQTL_start", 
				z = 'Z_score') %>%
			select(-"...10")


z -> p -> 1e-6

2e-7


/sc/arion/projects/CommonMind/zengb02/single_cell_eQTL_for_Gabriel_02212024/eQTL_detection_on_pearsonclass_level/trans-eQTL_detection/another/for_MAPT_region/analysis/significant_trans_eQTL_in_*

fig = df_top %>%
	select(Cell, Gene) %>%
	group_by(Cell) %>%
	summarize(nGenes = length(unique(Gene))) %>%
	mutate(Cell = factor(Cell, ord)) %>%
	ggplot(aes(nGenes, Cell, fill=Cell, label=nGenes)) +
		geom_bar(stat="identity") +
		scale_x_continuous(limits=c(0, 400), expand=c(0,0)) +
		theme_classic() +
		theme(aspect.ratio=1, legend.position="none") +
		scale_fill_manual(values = cols) +
		geom_text()

ggsave(fig, file="plots/trans_eQTL_count.pdf")

res = lapply( unique(df_top$Cell), function(CT){
	df_top %>%
		filter(Cell == CT) %>%
		pull(Gene) %>%
		unique
	})
names(res) = unique(df_top$Cell)


ord = ord[ord %in% names(res)]
pdf("plots/trans_upset.pdf", height=5, width=5)
upset(fromList(res[ord]), order.by = "freq", nsets=8, keep.order=TRUE, sets=ord, nintersects=15)
dev.off()

df_top %>%
	filter(Gene == "CBWD1")




df_top %>%
	arrange(-abs(z))

##############
# MAPT locus #
##############


cd /sc/arion/projects/CommonMind/hoffman/sceqtl/sumStats/
R


library(ggplot2)
library(tidyverse)
library(arrow)

files = system("ls /sc/arion/projects/CommonMind/zengb02/single_cell_eQTL_for_Gabriel_02212024/eQTL_detection_on_pearsonclass_level/trans-eQTL_detection/another/for_MAPT_region/analysis/significant_trans_eQTL_in_*", intern=TRUE)



df = lapply(files, function(file){
	dfin = read_tsv(file, col_names = FALSE, show_col_types=FALSE, progress=FALSE, num_threads=1)
	colnames(dfin) = c("chrom", 'start', 'end', 'variant', 'allele1', 'allele2_assessed', 'gene', 'beta', 'se', 'z')

	# CellType = gsub("eQTL_result_MSSM_(\\S+)_.*_.*$", "\\1", basename(file))
	# dfin$CellType = CellType
	# dfin
})








# Oct 21
#######

# cd /Users/gabrielhoffman/Library/CloudStorage/Dropbox/projects/PsychAD3.0/supplementary_figures/trans-eQTL/circos_plot

library(ggplot2)
library(tidyverse)

# hubs
df = read.table("adjust_trans_eQTL_hits_hub") 
colnames(df) = c("CellType", "SNP", "N")

cols = c("Immune" = "#c70813", "Astro" = "#d2af81", "OPC" = "#8f8a80","Oligo"="#d5e4a2", "EN"="#197ec0","IN"="#1a9993")

fig = df %>%
	group_by(CellType, N) %>%
	summarize(Count = n()) %>%
	mutate(CellType = factor(CellType, rev(names(cols)))) %>%
	ggplot(aes(factor(N), Count, fill=CellType, label=Count)) +
		geom_bar(stat="identity") +
		facet_wrap(~CellType, nrow=1) +
		theme_classic() +
		theme(aspect.ratio=1, legend.position="none") +
		scale_fill_manual("Cell", values = cols) +
		scale_y_continuous(limits=c(0,24), expand=c(0,0)) +
		geom_text() +
		xlab("Number of trans genes in regulatory hub")
 
ggsave(fig, file="trans_hubs.pdf")


# Trans/coloc overlap
#####################

# cd /Users/gabrielhoffman/Library/CloudStorage/Dropbox/projects/PsychAD3.0/supplementary_figures/trans-eQTL/Trans-eQTL_GWAS_integration

library(tidyverse)

df = read_tsv("merged_infor_for_trans_eGenes_6_cells")



df %>% 
  filter(ppH4 > 0.8) %>%
  filter(cell_trans_ == cell_eQTL)



df %>% 
  filter(ppH4 > 0.5) %>%
  pull(gene) %>%
  unique %>%
  length




ord = c("EN", "IN",  "Oligo", "OPC", "Astro", "Immune")



# heatmap sorted by genes
#########################

for(ds in names(table(df$disease)) ){
  df2 = df %>% 
  	# filter(gene %in% genes) %>%
    mutate(cell_eQTL = factor(cell_eQTL, ord)) %>%
    mutate(cell_trans_ = factor(cell_trans_, ord)) %>%
    filter(!is.na(cell_eQTL)) %>%
    complete(cell_trans_, cell_eQTL, gene, disease, fill=list(ppH4 = 0)) %>%
    filter(disease == ds) 

  # subset to genes with a non-zero score
  genes = df2 %>%
    group_by(gene) %>%
    summarize(max = max(ppH4)) %>%
    filter(max > 0) %>%
    pull(gene)

  if( length(genes) > 3){
    # order
    M = df2 %>%
      filter(gene %in% genes) %>%
      mutate(cell2 = paste(cell_trans_, cell_eQTL)) %>%
      select(-disease, -cell_trans_, -cell_eQTL) %>%
      pivot_wider(names_from = gene, values_from = ppH4 ) %>%
      column_to_rownames('cell2') %>%
      as.matrix %>%
      t

    hcl = hclust( dist(M))
    genes = hcl$labels[hcl$order]
  }

  df2$txt = ""
  df2$txt[df2$ppH4 > 0.5] = "*"
  df2$txt[df2$ppH4 > 0.8] = "#"

  fig = df2 %>%
    filter(gene %in% genes) %>%
    mutate(gene = factor(gene, genes)) %>%
    ggplot(aes(cell_eQTL, gene, color=ppH4, size=ppH4, label=txt)) +
      geom_point(stroke=0) +
      facet_grid(disease~cell_trans_) +
      coord_equal() +
      scale_color_gradient(low="white", high="red", limits=c(0,1)) +
      scale_size_area(limits=c(0,1)) +
      theme_classic() +
      scale_x_discrete(guide = guide_axis(angle = 90), name="Cell class with colocalization signal") +
      ggtitle(ds) +
      geom_text( aes( vjust=0.5), color="black")

  file = paste0("coloc_trans_overlap_", ds, ".pdf")
  ggsave(fig, file=file, height=7, width=9.5)
}


  










