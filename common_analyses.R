



plot_within_btw = function( simMat, design, col=c("gray", "steelblue1"),main='',ylim=NULL,...){
  variation = list(type=c(), value=c(), donor1=c(), donor2=c())
  for(i in 1:ncol(design)){
    idx1 = which(design[,i] == 1)   
    idx0 = which(design[,i] == 0)
    C_sub = simMat[idx1, idx1]
    variation$type = append( variation$type, rep("Same Donor", sum(lower.tri(C_sub))) )  
    variation$value = append( variation$value, C_sub[lower.tri(C_sub)])   
    variation$donor1 = append(variation$donor1, rep(colnames(design)[i], length(C_sub[lower.tri(C_sub)]) ))      
    variation$donor2 = append(variation$donor2, rep(colnames(design)[i], length(C_sub[lower.tri(C_sub)]) ))
    variation$type = append( variation$type, rep("Different Donor", length(simMat[idx1,idx0])) )  
    variation$value = append( variation$value, simMat[idx1,idx0])    
    variation$donor1 = append(variation$donor1, rep("", length(simMat[idx1,idx0])) )      
    variation$donor2 = append(variation$donor2, rep("", length(simMat[idx1,idx0])))
  }
  variation = as.data.frame(variation)
  fig = ggplot( variation, aes(type, value)) + 
    geom_violin(scale = "width", aes(fill=type)) + 
    ylab("Variance explained (%)") + 
    xlab("") + 
    geom_boxplot(width = 0.07, fill = "grey", outlier.color = "black") + 
    theme_bw(15) + 
    scale_fill_manual(values = col) + 
    theme(legend.position = "none") + 
    ylab("Correlation between samples") 

    if( main !='' ){
      fig = fig + ggtitle( main )
    }
    if( !is.null(ylim) ){
      fig = fig + ylim( ylim )
    }
    fig = fig + theme(aspect.ratio=1)
    return(list(fig=fig, variation=variation))
}


eval_within_across_donor = function( pbObj, Donor ){

  formula = as.formula(paste0("~ 0 + ", Donor))

	# for each assay, evaluate within/across donor correlation
	resList = lapply( assayNames(pbObj), function(k){

	  # extract gene exrpression data
	  geneExpr = assay(pbObj, k)

	  # Evaluate correlation
	  C = cor(geneExpr$E, method="spearman")

	  # extract data
	  info = colData(pbObj)[rownames(C),,drop=FALSE]

	  # Evalute within-donor similarity
	  res = plot_within_btw(C, model.matrix(formula, info))

	  # return correlations
	  data.frame(assay = k, res$variation)
	})
	df = do.call(rbind, resList)

	# combind plots
	ggplot( df, aes(type, value)) + geom_violin(scale = "width", aes(fill=type)) + 
	    ylab("Variance explained (%)") + xlab("") + geom_boxplot(width = 0.07, 
	        fill = "grey", outlier.color = "black") + theme_bw(15) + 
	        scale_fill_manual(values = col) + theme(legend.position = "none", aspect.ratio=1, axis.text.x = element_text(angle = 45, hjust=1, vjust=1)) + 
	        ylab("Correlation between samples") + scale_fill_manual(values=c("grey", "steelblue1")) + facet_wrap(~assay)
}



# Compute correlation and standard error for Pearson and Spearman
cor.se = function(x,y, method = c("pearson", "kendall", "spearman"),...){

  method = match.arg(method)

  if( method == "pearson"){
    df <- cor.test(x,y, method=method,...) %>%
      tidy %>%
      mutate(se = sqrt((1 - estimate^2)/parameter))
    res = data.frame(rho = df$estimate, rho.se = df$se)
  }else if(method == "spearman"){
    # https://stats.stackexchange.com/questions/18887/how-to-calculate-a-confidence-interval-for-spearmans-rank-correlation
    rho <- cor(x,y, method=method,...)
    n <- sum(complete.cases(x, y))
    rho.se <- sqrt((1+rho^2/2)/(n-3))

    res <- data.frame(rho = rho, rho.se = rho.se)
  }

  return(res) 
}

