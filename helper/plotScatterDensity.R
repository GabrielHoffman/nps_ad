
library(ggplot2)
library(viridis)

# Compute two dimensional density
get_density <- function(x, y, n = 250) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

# make scatter plot with 2D density using viridis colors
plotScatterDensity = function(value1, value2){

  # convert two vectors to a data frame
	df = data.frame(cbind(value1, value2))

  # determine limits of the plot
	lim = with(df, max(abs(c(value1, value2))))
  
  # Compute 2D density
	df$density <- get_density(df$value1, df$value2, n = 100)
	
  # Scatter plot colored by density
  ggplot(df, aes(value1, value2, color=density)) + geom_point(size=.4) + theme_bw(16) + theme(aspect.ratio=1, legend.position="bottom", plot.title = element_text(hjust = 0.5)) + geom_abline(color="red") + geom_vline(xintercept=0, col="grey40", linetype="dashed") + geom_hline(yintercept=0, col="grey40", linetype="dashed") + xlim(-lim, lim) + ylim(-lim, lim) + scale_color_viridis() + geom_smooth(method="lm", se=FALSE, color="darkorange")
}
