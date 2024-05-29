


# write identifers for each cell

library(SingleCellExperiment)
library(zellkonverter)
library(tidyverse)

path = "/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/h5ad_final/"

h5ad_files = c(
RUSH = "RUSH_2024-02-01_14_53.h5ad", 
HBCC = "HBCC_2024-02-01_15_18.h5ad",
MSSM = "MSSM_2024-02-01_16_17.h5ad")
h5ad_files = paste0(path, h5ad_files)

df = lapply(h5ad_files, function(file){
	sce = readH5AD(file, use_hdf5=TRUE, verbose=FALSE, version="0.8.0")

	colData(sce) %>%
			as.data.frame %>%
			rownames_to_column %>%
			as_tibble %>%
			mutate(CellID = rowname) %>%
			select(CellID, class, subclass, subtype, Source)
	})
df = bind_rows(df)


file = paste0(path, "cell_identifiers.tsv.gz")
write_tsv(df, file=file)



