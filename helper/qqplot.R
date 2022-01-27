
library(ggplot2)

qqplot = function( pvals, pmin=1e-300){

	N = length(pvals)
	ci = 0.95
	df = data.frame(x 		= ppoints(N),
					pvals 	= pmax(pmin, sort(pvals)),
				    lower   = qbeta(p = (1 - ci) / 2, shape1 = 1:N, shape2 = N:1, lower.tail=FALSE),
				    upper   = qbeta(p = (1 + ci) / 2, shape1 = 1:N, shape2 = N:1, lower.tail=FALSE))

	fig = ggplot(df, aes(-log10(x), -log10(pvals))) + geom_ribbon(aes(ymin = -log10(lower), ymax = -log10(upper)), fill = "grey90") 

	xmax = -log10(min(df$x))
	ymax = -log10(min(c(df$upper, df$pvals, 0.05/N * 0.5)))

	a = data.frame(x=c(0,xmax), y=c(0, xmax))

	fig = fig + geom_line(aes(x,y), color="red", data=a) + scale_x_continuous(limits=c(0, xmax*1.02), expand=c(0,0)) + scale_y_continuous(limits=c(0, ymax*1.02), expand=c(0,0))

	fig + geom_point(size=.5) + xlab(bquote(Expected~-log[10]~P)) + ylab(bquote(Observed~-log[10]~P)) + theme_classic(12) + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5))
}

