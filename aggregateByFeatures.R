
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

  if( length(chrom) !+ length(feature) ){
    stop("chrom and feature must be the same length")
  }

  if( ! assay %in% assayNames(sceCombine) ){
    stop("assay not found in SingleCellExperiment: ", assay)
  }

  # keep only chroms/features that are in sce
  idx = feature %in% rownames(sce)
  feature = feature[idx]
  chroms = chrom[idx]

  # for each chromosome
  chromExpr = bplapply( as.character(sort(unique(chroms))), function(chrom, counts){

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
    idx = unique(c(seq(1, ncol(counts), by=batchSize), ncol(counts)))

    chromCounts = lapply( 2:length(idx), function(i){
      colSums2( counts[feature_local,idx[i-1]:idx[i],drop=FALSE] )
    })
    chromCounts = do.call(cbind, chromCounts)

    browser()
    
    # return as data.frame
    df = data.frame(chromCounts, check.names=FALSE)
    rownames(df) = chrom
    colnames(df) = colnames(sce)

    df
    }, counts = assay(sce, assay), BPPARAM=BPPARAM)
  chromExpr = do.call(rbind, chromExpr)

  chromExpr
}

# if( ! chromCol %in% colnames(geneInfo) ){
#     stop("geneInfo must have column ", chromCol)
#   }
#   if( ! geneCol %in% colnames(geneInfo) ){
#     stop("geneInfo must have column ", geneCol)
#   }