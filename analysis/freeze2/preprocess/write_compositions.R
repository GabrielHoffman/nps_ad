# Gabriel Hoffman
# June 28, 2023
# Write cell compositions for each donor to file

library(dreamlet)
library(crumblr)

dataset = c("RUSH", "MSSM", "AGING", "HBCC", "FULL")
cluster_id_options = c("class", "subclass", "subtype")

df = expand.grid(dataset = dataset, cluster_id = cluster_id_options)

res = lapply(seq(nrow(df)), function(i){

	cat("\r", i, '    ')
	# Read pseudobulk
	folder = "/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/pseudobulk/"
	pattern = paste0("", df$dataset[i], "_2024-02-01_.*_PB_SubID_", df$cluster_id[i], ".RDS")
	file = dir(folder, pattern=pattern, full.names=TRUE)

	if( length(file) == 0) stop("No file found")
	if( length(file) > 2) stop("Too many files")
	if( length(file) == 2){			
		a = file.info(file)$ctime
		file = ifelse(difftime(a[1], a[2]) > 0, file[1], file[2])
	}

    pb = readRDS(file)

    # crumblr transform
    cobj = crumblr(cellCounts(pb))

    # CLR
    outfile = paste0(dirname(folder), "/composition/", gsub("_PB_", "_clr_", gsub("RDS$", "tsv", basename(file))))
    write.table(cobj$E, file=outfile, quote=FALSE, sep="\t")
    R.utils::gzip(outfile, overwrite=TRUE)

    # Weighted CLR
    fit = dream(cobj, ~ 1, colData(pb))
    resid = residuals(fit, cobj, type="pearson")

     outfile = paste0(dirname(folder), "/composition/", gsub("_PB_", "_clrPearson_", gsub("RDS$", "tsv", basename(file))))
    write.table(resid, file=outfile, quote=FALSE, sep="\t")
    R.utils::gzip(outfile, overwrite=TRUE)
	})


