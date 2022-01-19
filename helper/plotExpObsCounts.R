
library(ggplot2)
library(viridis)

plotExpObsCounts = function(ods, geneID, log=TRUE){

  # Compute two dimensional density
  get_density <- function(x, y, n = 250) {
    dens <- MASS::kde2d(x = x, y = y, n = n)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
    return(dens$z[ii])
  }

  ods <- ods[geneID]
  df <- data.table(
          feature_id = geneID,
          sampleID   = colnames(ods),
          observed   = as.vector(counts(ods)) + isTRUE(log),
          expected   = as.vector(normalizationFactors(ods)) + isTRUE(log),
          aberrant   = as.vector(aberrant(ods)))

  # determine limits of the plot
  lim = max(max(df$expected), max(df$observed))

  # Compute 2D density
  df_null = df[!df$aberrant,]
  df_null$density <- with(df_null, get_density(log10(expected), log10(observed), n = 160))

  # Scatter plot colored by density
  fig = ggplot(df_null, aes(log10(expected), log10(observed), color=density)) + geom_point(size=.03) + theme_classic(16) + theme(aspect.ratio=1, legend.position="bottom", plot.title = element_text(hjust = 0.5)) + scale_x_continuous(limits=c(1, log10(lim)), breaks=0:5, labels=10^c(0:5)) + scale_y_continuous(limits=c(1, log10(lim)), breaks=0:5, labels=10^c(0:5) ) + scale_color_viridis() + xlab(paste('Expected counts', ifelse(isTRUE(log), '+ 1', ''))) + ylab(paste('Observed counts', ifelse(isTRUE(log), '+ 1', ''))) + ggtitle(geneID)

  down = with(df, which(aberrant & observed < expected))
  fig2 = fig + geom_point(aes(log10(expected), log10(observed)), data=df[down,], color="dodgerblue", size=.3, stroke=.3)

  up = with(df, which(aberrant & observed > expected))
  fig3 = fig2 + geom_point(aes(log10(expected), log10(observed)), data=df[up,], color="#ff1e20", size=.3, stroke=.3)

  fig4 = fig3 + annotate("text", x = c(1.12), y=c(log10(lim)*.9), label=c(paste("Up:", length(up), "\nDown:", length(down))), hjust = 0)  

  fig4
}

