

library(dreamlet)
library(SingleCellExperiment)
library(tidyverse)
library(zenith)
library(qvalue)
library(cowplot)
library(synapser)
synLogin()

# load(synGet("syn51114763")$path)
load("/hpc/users/hoffmg01/contrasts_for_dreamlet.Rdata")

res.proc = readRDS("/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/processAssays/HBCC_2023-02-28_12_36_processAssays_Channel_subclass.RDS")



# Set base covariates for dream; do not add "Brain_bank" (if samples are from multiple brain banks, this covariate will be included in $covariates_incl)
covariates_base = c("Age", "Sex", "PMI", "(1|SubID)", "(1|poolID)")

# TODO / FIXME: make sure that the loaded metadata (data frame `metadata`) replaces/updates the metadata attached to SingleCellExperiment - as those don't have contrasts-specific target columns such as *c01x* etc.

# Run the analysis for categorical contrasts
for( ctr in  names(CONTRASTS$CAT[[8]] ){

	ctr = CONTRASTS$CAT[[8]]
	print(ctr$name)
	print(paste0("> Variable: ", ctr$variable))
	print(paste0("> Vector of contrasts: ", paste0(ctr$contrasts, collapse=",")))

	# Set covariates: from base covars, add extra covars relevant for given contrast (e.g. "Brain_bank" for contrasts mixing samples from different cohorts) and remove covars that could collide with variable (e.g. Age for "Aging" contrast)
	covariates = setdiff(union(covariates_base, ctr$covariates_incl), ctr$covariates_excl)

	# variable you are doing the test on
	variable = ctr$variable

	# name results and define contrasts
	contrasts = ctr$contrasts

	# user-specified formula
	form = paste0("~ 0 + ", variable, " + ", paste0(covariates, collapse=" + "))
	form = as.formula(form)

	# merge subject-level metadata with colData(res.proc)
	metadata_sub = metadata[metadata$SubID %in% colData(res.proc)$SubID,]
	idx = match( colData(res.proc)$SubID, metadata_sub$SubID)

	# df = cbind(as.character(colData(res.proc)$SubID), 
	# 	metadata_sub$SubID[idx])
	# all.equal(df[,1], df[,2])

	cols = colnames(metadata_sub) %in% colnames(colData(res.proc))
	# colData(res.proc) = cbind(colData(res.proc), metadata_sub[idx,!cols])
	res.proc@data = cbind(colData(res.proc), metadata_sub[idx,!cols])

	# run dreamlet
	fit = dreamlet( res.proc, form, contrasts = contrasts, BPPARAM=SnowParam(4, progressbar=TRUE))

	saveRDS(fit, file="~/www/BD/fit.RDS")
}

# cell type order
ctorder = c('EN_L2_3_IT', 'EN_L3_5_IT_1', 'EN_L3_5_IT_2', 'EN_L3_5_IT_3', 'EN_L5_6_NP', 'EN_L6_CT', 'EN_L6_IT', 'EN_NF', 'IN_ADARB2', 'IN_LAMP5', 'IN_PVALB', 'IN_PVALB_CHC', 'IN_SST', 'IN_VIP', 'Oligo', 'OPC', 'Astro', 'Micro_PVM', 'CD8_T', 'PC', 'VLMC','Endo', "Immune")

ctorder = ctorder[ctorder %in% assayNames(fit)]

coef = "c08xBD - c08xControl"


path = "~/www/BD/"

file = paste0(path, "plotVolcano", ".pdf")
fig = plotVolcano( fit, coef = coef, assays=ctorder )
ggsave(file, fig, height=20)


tab = topTable(fit, coef = coef, number=Inf )

res = tab %>%
	as_tibble %>%
	filter(adj.P.Val < 0.05) %>%
	group_by(ID) %>%
	summarise(nDE = length(ID)) %>%
	filter(nDE > 2) %>%
	arrange(desc(nDE), ID)


genes = rev(sort(res$ID))
file = paste0(path, "plotGeneHeatmap", ".pdf")
fig = plotGeneHeatmap( fit, coef=coef, genes=genes, assays=ctorder)
ggsave(file, fig)


df_de = tab %>%
	as_tibble %>%
	group_by(assay) %>%
	summarise( 
	  nGenes = length(adj.P.Val), 
	  nDE = sum(adj.P.Val < 0.05),
	  pi1 = 1 - qvalue(P.Value)$pi0) %>%
	mutate(assay = factor(assay, ctorder)) 


ymax = 1.05*max(df_de$nGenes)
fig1 = ggplot(df_de, aes(nGenes, assay, fill=assay)) + 
    geom_bar(stat="identity") + 
    theme_classic() +
    theme(aspect.ratio=1, legend.position="none") +
    scale_x_continuous(limits=c(0,ymax), expand=c(0,0)) +
    xlab("# genes expressed") +
    ylab("Cell type") 

ymax = max(1.05*max(df_de$nDE), 100)
fig2 = ggplot(df_de, aes(nDE, assay, fill=assay)) + 
    geom_bar(stat="identity") + 
    theme_classic() +
    theme(aspect.ratio=1, legend.position="none", axis.text.y=element_blank()) +
    scale_x_continuous(limits=c(0,ymax), expand=c(0,0)) +
    xlab("# genes with FDR < 5%") +
    ylab('')

fig3 = ggplot(df_de, aes(pi1, assay, fill=assay)) + 
    geom_bar(stat="identity") + 
    theme_classic() +
    theme(aspect.ratio=1, legend.position="none", axis.text.y=element_blank()) +
    scale_x_continuous(limits=c(0,1), expand=c(0,0)) +
    xlab(bquote(pi[1]))+
    ylab('')

fig = plot_grid(fig1, fig2, fig3, labels=LETTERS[1:3], nrow=1, axis="tblr", align="hv")

file = paste0(path, "nDE", ".pdf")
ggsave(file, fig, width=12)

figList = lapply(sort(res$ID), function(gene){
	plotForest( fit, coef = coef, gene = gene) + theme(legend.position="right", aspect.ratio=1)
})
fig = plot_grid(plotlist=figList, ncol=3)

file = paste0(path, "plotForest", ".pdf")
ggsave(file, fig, height=20, width=12)



go.gs = get_GeneOntology(to="SYMBOL")
   
# Run zenith gene set analysis on result of dreamlet
res_zenith = zenith_gsa(fit, coef = coef, go.gs)


fig = plotZenithResults(res_zenith, 5, 1)


file = paste0(path, "zenith", ".pdf")
ggsave(file, fig, height=20, width=12)



res_mash = run_mash(fit, coef=coef)



df_gs = zenith_gsa(res_mash, go.gs)

# Heatmap of results
fig = plotZenithResults(df_gs, 5, 1)


file = paste0(path, "zenith_mashr", ".pdf")
ggsave(file, fig, height=20, width=12)









