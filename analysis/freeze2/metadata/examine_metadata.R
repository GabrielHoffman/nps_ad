

BRAAK_AD
Plq_Mn
CERAD
CDRScore
Plq_Mn_MFG
Cognitive_Resilience

cn = colnames(colData(pb))[29:68]
# cn = cn[ !cn %in% 'Vascular']
info = data.frame(colData(pb)[,cn,drop=FALSE])
info[is.na(info)] = 0

ph = apply(info, 1, function(x) paste(x[x!="no"], collapse="/"))

sort(table(ph))


library(dreamlet)

path = "/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/h5ad_final/"

h5ad_files = c(
RUSH = "RUSH_2023-09-12_16_13.h5ad", 
HBCC = "HBCC_2023-09-12_16_28.h5ad",
MSSM = "MSSM_2023-09-12_17_04.h5ad",
FULL = "FULL_2023-09-12_18_32.h5ad",
AGING = "AGING_2023-09-12_20_35.h5ad")

h5ad_files = sapply(h5ad_files, function(x) paste0(path, x))


sample_id = "SubID"
cluster_id = "class"
pbList = lapply( names(h5ad_files[-4]), function(id){

	params = list(DATASET = h5ad_files[id])

	filePrefix = paste0("/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/pseudobulk/", gsub(".h5ad$", "", basename(params$DATASET)), "_PB_", sample_id, "_", cluster_id)
	file = paste0(filePrefix, ".RDS")
	readRDS(file)
})
names(pbList) = names(h5ad_files[-4])



augment_data = function(pb){

	cn1 = c( 'MCI', 'SCZ', 'AD')
	cn2 = c( 'DLBD', 'Dementia', 'FTD')
	cn_bd = c("BD_unspecific", "BD_I", "BD_II")
	cn_all = c(cn1, cn2, cn_bd)

	info = data.frame(colData(pb)[,cn_all])
	info[is.na(info)] = 0

	# code BD
	info$BD = apply(info[,cn_bd], 1, function(x) ifelse(length(unique(x)) == 1, 0, 1))

	for(id in colnames(info)){
		info[,id] = factor(c("no", id)[info[,id]+1], levels=c("no", id))
	}

	# set 1
	cn1 = c(cn1, "BD")
	ph1 = apply(info[,cn1], 1, function(x) paste(x[x!="no"], collapse="/"))
	ph1[ph1 == "SCZ/AD"] = "AD"
	ph1[ph1 == "MCI/SCZ"] = "SCZ"
	ph1[ph1 == ""] = "none"
	# sort(table(ph1))
	pb$MajorDisease = factor(ph1)

	# set2
	ph2 = apply(info[,cn2], 1, function(x) paste(x[x!="no"], collapse="/"))
	ph2[ph2 == "Dementia/FTD"] = "FTD"
	ph2[ph2 == "DLBD/Dementia"] = "DLBD"
	ph2[ph2 == ""] = "none"
	# sort(table(ph2))
	pb$DementiaType = factor(ph2)

	if( length(unique(pb$DementiaType)) == 1){
		pb$DementiaType = c()
	}
	if( length(unique(pb$MajorDisease)) == 1){
		pb$MajorDisease = c()
	}

	# metrics
	if( any(!is.na(pb$BRAAK_AD)) ){
		pb$mod_BRAAK_AD = factor(as.character(pb$BRAAK_AD))
		pb$mod_BRAAK_AD[is.na(pb$mod_BRAAK_AD)] = 2
	}

	if( any(!is.na(pb$Plq_Mn)) ){
		pb$mod_Plq_Mn = pb$Plq_Mn
		pb$mod_Plq_Mn[is.na(pb$mod_Plq_Mn)] = mean(pb$mod_Plq_Mn, na.rm=TRUE)
	}

	if( any(!is.na(pb$CERAD)) ){	
		pb$mod_CERAD = pb$CERAD
		pb$mod_CERAD[is.na(pb$mod_CERAD)] = mean(pb$mod_CERAD, na.rm=TRUE)
	}

	if( any(!is.na(pb$CDRScore)) ){
		pb$mod_CDRScore = pb$CDRScore
		pb$mod_CDRScore[is.na(pb$mod_CDRScore)] = mean(pb$mod_CDRScore, na.rm=TRUE)
	}

	if( any(!is.na(pb$Cognitive_Resilience)) ){
		pb$mod_Cognitive_Resilience = pb$Cognitive_Resilience
		pb$mod_Cognitive_Resilience[is.na(pb$mod_Cognitive_Resilience)] = mean(pb$mod_Cognitive_Resilience, na.rm=TRUE)
	}

	pb
}

a = augment_data(pb)

form = ~ Sex + scale(Age) + log(n_genes) + (1|mod_BRAAK_AD) + mod_Plq_Mn + (1|mod_CERAD) + (1|mod_CDRScore) + mod_Cognitive_Resilience + (1|MajorDisease) + (1|DementiaType) + (1|cogdx)


# Sex, Age, log(n_genes)
# Braak: categorical
# Plaque mean: continuous
# CERAD: categorical
# CRD: categorical in MSSM
# cogdx: categorical in RUSH
# Cognitive_Resilience: continuous 
# MajorDisease: AD, MCI, SCZ, BD, none
# DementiaType: Dementia, DLBD, FTD, none

res.proc = processAssays( pb[1:100,], form, BPPARAM = SnowParam(4), assays="Astro")




form = ~ Sex + scale(Age) + scale(PMI) + (1|mod_BRAAK_AD) + mod_Plq_Mn + (1|mod_CERAD) + (1|mod_CDRScore) + mod_Cognitive_Resilience + (1|MajorDisease) + (1|DementiaType) + log(n_genes) 

pdf("~/www/test.pdf")
figList = lapply( names(pbList), function(id){

	message(id)

	pb.mod = augment_data(pbList[[id]])

	i = (all.vars(form) %in% colnames(colData(pb.mod))) 
	vars.keep = all.vars(form)[i]

	j = apply(colData(pb.mod)[,vars.keep,drop=FALSE], 2, function(x) length(unique(x)) > 1)
	vars.keep = vars.keep[j]

	form2 = paste("~", paste(vars.keep , collapse=" + "))

	C = canCorPairs(form2, colData(pb.mod))

	rownames(C) = gsub("mod_", "", rownames(C))
	colnames(C) = gsub("mod_", "", colnames(C))

	plotCorrMatrix(C, main=id)
})
dev.off()













