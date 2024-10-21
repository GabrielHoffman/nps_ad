

/sc/arion/projects/CommonMind/zengb02/single_cell_eQTL_for_Gabriel_02212024/plot/plot_for_cell_type_specific_eQTL/run_plot_for_caQTL_eQTL_GWAS_colocalization_v3.R

cd /hpc/users/hoffmg01/work/tmp

/sc/arion/projects/CommonMind/zengb02/single_cell_eQTL_for_Gabriel_02212024/plot/plot_for_cell_type_specific_eQTL/run_plot.sh test


R --args APP Peak_132559 Oligo neuron AD3




plotPosterior2 = function( grObj, wh, size=8 ){
  	gr_sub = grObj[grObj %within% wh]
	Max = 0
	if(length(score(gr_sub))>0){
	 Max<-max(score(gr_sub))+1
	}
	
	df = data.frame( gr_sub )

	ggplot(df, aes(start, score, size=pip)) + 
			geom_point(aes(color=pip)) +
			theme_bw(size) +
			scale_color_gradient2(low="black", mid="orange", high="red", limits=c(1e-4,1), na.value="black") + 
			ylab(bquote(-log[10]~P)) + 
			scale_y_continuous(expand=c(0,0.01), limits=c(0, Max)) + 
			scale_x_continuous(expand=c(0,0), limits=c(start(wh), end(wh)), label=comma_format()) + 
			theme(legend.position="none",axis.line.x=element_line(size=.5, colour="black"),
			 axis.line.y=element_line(size=.5, colour="black"),
			 axis.ticks=element_line(color="black"),axis.text=element_text(color="black", size=7),
			 axis.title=element_text(size=7)) 
}



fig_mht_eQTL1 = plotPosterior2(gr_eQTL11_sub, wh)

ggplot2::ggsave(fig_mht_eQTL1, file="test.png")




df_sub[which.max(df_sub$pip),]


gr_eQTL11_sub<-gr_eQTL1_sub

idx = match(gr_eQTL11_sub$Variant, gr_eQTL1_finemap_sub$Variant)

i1 = which(!is.na(idx))
i2 = idx[which(!is.na(idx))]

gr_eQTL11_sub$pip = 0
gr_eQTL11_sub$pip[i1] = gr_eQTL1_finemap_sub$score[i2]
 
plotPosterior2 = function( grObj, wh, size=8 ){
  	gr_sub = grObj[grObj %within% wh]
	Max = 0
	if(length(score(gr_sub))>0){
	 Max<-max(score(gr_sub))+1
	}
	
	df = data.frame( grObj )

	ggplot(df, aes(start, score, size=pip)) + 
			geom_point(aes(color=pip)) +
			theme_bw(size) +
			scale_color_gradient2(low="black", mid="orange", high="red", limits=c(1e-4,1), na.value="black") + 
			ylab(bquote(-log[10]~P)) + 
			scale_y_continuous(expand=c(0,0.01), limits=c(0, Max)) + 
			scale_x_continuous(expand=c(0,0), limits=c(start(wh), end(wh)), label=comma_format()) + 
			theme(legend.position="none",axis.line.x=element_line(size=.5, colour="black"),
			 axis.line.y=element_line(size=.5, colour="black"),
			 axis.ticks=element_line(color="black"),axis.text=element_text(color="black", size=7),
			 axis.title=element_text(size=7)) 
}



