
#' Sum counts across chromosome
#'
#' Sum counts across chromosome
#'
#' @param sce SingleCellExperiment
#' @param assay which assay to summarize
#' @param geneInfo \code{data.frame} storing gene ids and chromosome
#' @param chromCol name of column in \code{geneInfo} storing chromosome
#' @param geneCol name of column in \code{geneInfo} storing gene
#' @param batchSize number of columns in a batch
#'
#' @import SingleCellExperiment
#' @importFrom DelayedMatrixStats colSums2 
aggregateByFeatures = function(sce, assay, chroms, feature, BPPARAM=SerialParam(), batchSize=100000){

  if( length(chroms) != length(feature) ){
    stop("chroms and feature must be the same length")
  }

  if( ! assay %in% assayNames(sceCombine) ){
    stop("assay not found in SingleCellExperiment: ", assay)
  }

  # keep only chroms/features that are in sce
  idx = feature %in% rownames(sce)
  feature = feature[idx]
  chroms = chroms[idx]

  chrUnique = as.character(sort(unique(chroms)))

  # divide columns into chunks
  d = seq(1, ncol(sce))
  idx = split(d, ceiling(seq_along(d)/batchSize))

  # for each chromosome
  chromExpr = bplapply( chrUnique, function(chrom, counts){

    suppressPackageStartupMessages({
    library(DelayedArray)
    library(HDF5Array)
    library(DelayedMatrixStats)
    })

    # features on specified chromsome
    feature_local = feature[chroms == chrom]

    message("\r", chrom, ": ", length(feature_local), '     ')

    # sum reads across genes on this chromosome
    # colSums2() causes overflow when used directly
    # here is a workaround
    # use colSums2() in batches
    chromCounts = lapply( idx, function(i){
      colSums2( counts[feature_local,i,drop=FALSE] )
    })
    chromCounts = do.call(c, chromCounts)

    chromCounts
    }, counts = assay(sce, assay), BPPARAM=BPPARAM)

  chromExpr = do.call(rbind, chromExpr)
  rownames(chromExpr) = chrUnique
  colnames(chromExpr) = colnames(sce)

  chromExpr
}

