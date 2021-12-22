
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
aggregateByFeatures = function(sce, assay, geneInfo, chromCol, geneCol, BPPARAM=SerialParam(), batchSize=100000){

  if( ! chromCol %in% colnames(geneInfo) ){
    stop("geneInfo must have column ", chromCol)
  }
  if( ! geneCol %in% colnames(geneInfo) ){
    stop("geneInfo must have column ", geneCol)
  }
   if( ! assay %in% assayNames(sceCombine) ){
    stop("assay not found in SingleCellExperiment: ", assay)
  }

  # unique set of chroms
  chroms = as.character(sort(unique(geneInfo[[chromCol]])))

  # for each chromosome
  chromExpr = bplapply( chroms, function(chrom, counts){

    suppressPackageStartupMessages({
    library(DelayedArray)
    library(HDF5Array)
    library(DelayedMatrixStats)
    })

    # genes on specified chromsome
    genes = geneInfo[[geneCol]][geneInfo[[chromCol]] == chrom]

    message("\r", chrom, ": ", length(genes), '     ')

    # sum reads across genes on this chromosome
    # colSums2() causes overflow when used directly
    # here is a workaround
    # use colSums2() in batches
    idx = c(seq(1, ncol(counts), by=batchSize), ncol(counts))

    chromCounts = lapply( 2:length(idx), function(i){
      colSums2( counts[genes,idx[i-1]:idx[i],drop=FALSE] )
    })
    chromCounts = do.call(c, chromCounts)

    # return as data.frame
    df = data.frame(t(chromCounts), check.names=FALSE)
    rownames(df) = chrom
    colnames(df) = colnames(sce)

    df
    }, counts = assay(sce, assay), BPPARAM=BPPARAM)
  chromExpr = do.call(rbind, chromExpr)

  chromExpr
}

