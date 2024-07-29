
cd /sc/arion/projects/CommonMind/zengb02/single_cell_eQTL_for_Gabriel_02212024/eQTL_detection_on_pearsonclass_level/cis-eQTL_detection/analysis/enrichment_of_QTL_in_epimarker/


library(ggplot2)
library(tidyverse)
library(qvalue)
library(scales)
#library(liftOver)
library(rtracklayer)



file = "/sc/arion/projects/psychAD/NPS-AD/freeze2_proc/230816_PsychAD_capstone_F1/230921_PsychAD_color_palette.csv"
df_colors = read_csv(file) %>%
				select(category, name, color_hex) %>%
				filter(category %in% c("class", "subclass", "subtype")) %>%
				mutate(Dataset = paste0(category, " / ", name)) %>%
				select(Dataset, color = color_hex) %>%
				bind_rows(tibble(Dataset = "bulk / bulk", color = "grey"))
cols = df_colors$color
names(cols) = df_colors$Dataset
names(cols) = gsub(".* / ", "", names(cols))


main_theme1 = theme_classic()+theme(panel.background=element_blank(),
                                    panel.grid=element_blank(),
                                    axis.line.x=element_line(size=.5, colour="black"),
                                    axis.line.y=element_line(size=.5, colour="black"),
                                    axis.ticks=element_line(color="black"),
                                    axis.text=element_text(color="black", size=7),
                                    axis.title=element_text(size=7),
                                    legend.position='none',
                                    legend.background=element_blank(),
                                    legend.key=element_blank(),
                                    legend.key.size = unit(1.2, "lines"),
                                    legend.text= element_text(size=7),
                                    text=element_text(family="Helvetica", size=7),aspect.ratio=1, plot.title = element_text(hjust = 0.5))

fdensity_files <- list.files('./results/corces/',full.names = TRUE)


fdensity_res <- tibble(fdensity_files) %>% mutate(file_content=map(fdensity_files,read_delim,delim=' ',col_names=FALSE)) %>% 
  unnest(file_content) %>% 
  mutate(name=basename(fdensity_files)) %>% 
  dplyr::select(-fdensity_files) %>% 
  separate(name,into=c('tmp','sig','n1','Annot'),sep='_') %>% 
  separate(Annot,into=c('annot','n2','n3','n4'),sep="\\.") %>% 
  dplyr::select(-tmp) %>% dplyr::select(-n1) %>% dplyr::select(-n2) %>% dplyr::select(-n3) %>% dplyr::select(-n4) %>% 
  mutate(sig=gsub('.sig5FDR.bed','',sig)) %>% 
  mutate(annot=gsub('.bed','',annot)) %>% 
  group_by(sig,annot) %>% 
  mutate(X3_norm=X3/sum(X3)) %>% 
  mutate(sig=gsub('\\.',' ',sig) %>% gsub('OPCs   COPs','OPCs / COPs',.)) %>% 
  mutate(annot=gsub('\\.',' ',annot)) %>% 
  mutate(annot=gsub(' corces','',annot)) %>% 
  mutate(annot=gsub('ExcitatoryNeurons','Excitatory neurons',annot)) %>% 
  mutate(annot=gsub('InhibitoryNeurons','Inhibitory neurons',annot)) %>% 
  mutate(annot=gsub('OPCs','OPCs / COPs',annot)) %>% 
  mutate(annot=gsub('specific','ATAC (specific)',annot)) %>% 
  mutate(sig=paste0(sig,' eQTL')) %>% 
  mutate(sig=factor(case_when(
    sig=='astro eQTL' ~ 'Astro',
  #  sig=='Endothelial cells eQTL' ~ 'Endo eQTL',
    sig=='glu eQTL' ~ 'EN',
  sig=='gaba eQTL' ~ 'IN',
  #  sig=='glia eQTL' ~ 'glia eQTL',
  #  sig=='Excitatory neurons eQTL' ~ 'Ex eQTL',
  #  sig=='Inhibitory neurons eQTL' ~ 'Inh eQTL',
    sig=='microg eQTL' ~ 'Micro',
    sig=='oligo eQTL' ~ 'Oligo',
    sig=='OPC eQTL' ~ 'OPC'),
  #  sig=='Pericytes eQTL' ~ 'Pericytes eQTL'),
    levels=c('Astro','EN','IN','Micro','Oligo','OPC'))) 
 # filter(grepl('ATAC',annot)) %>% 
 # mutate(annot=factor(case_when(
 #   annot=='Astrocytes' ~ 'Astro ATAC',
    #annot=='Excitatory' ~ 'Ex ATAC',
    #annot=='Inhibitory' ~ 'Inh ATAC',
    #annot=='Microglia' ~ 'Micro ATAC',
    #annot=='Oligodendrocytes' ~ 'Oligo ATAC',
    #annot=='OPCs / COPs' ~ 'OPCs ATAC'),levels=c('Astro ATAC','Micro ATAC','Oligo ATAC','OPCs ATAC','Ex ATAC','Inh ATAC')))

p <- fdensity_res  %>% filter(!sig%in%c('Pericytes eQTL','Endo eQTL')) %>% ggplot(.,aes((X1+X2)/2,X3_norm,col=sig)) + geom_line() +
  facet_grid(sig~annot) + theme_classic() + theme(legend.position = 'none',text=element_text(size=18),axis.title = element_text(size=30)) + 
  xlab('Distance to eQTL (kb)') + ylab('Density') + scale_x_continuous(breaks=c(-500000,0,500000),labels = c('-500','0','500')) + 
  scale_y_continuous(breaks=c(0.005,0.010,0.015),labels = c('0.005','0.01','0.015')) + 
  #scale_color_manual(values=cols)+
  main_theme1+
  theme(text=element_text(size=10),axis.title = element_text(size=14))+
  scale_color_manual(breaks=c("Astro eQTL","EX-neuron eQTL","IN-neuron eQTL","Micro eQTL","Oligo eQTL","OPCs eQTL"),values=cols) +
  geom_vline(xintercept = 0,linetype="dashed",color="grey")
  
file = "/sc/arion/projects/CommonMind/hoffman/sceqtl/sumStats/plots/OCR_overlap.pdf"
ggsave(p,filename = file,width=8,height=11.5)







