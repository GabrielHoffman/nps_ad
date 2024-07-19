

cd /sc/arion/projects/CommonMind/hoffman/sceqtl/sumStats/
R 


library(tidyverse)
library(UpSetR)

files = system("ls /sc/arion/projects/CommonMind/zengb02/single_cell_eQTL_for_Gabriel_02212024/eQTL_detection_on_pearsonsubclass_level/cis-eQTL_detection/analysis/evaluate_eQTL_sharing_mashR/get_lfsr/lfsr_results_for_psychAD_sceQTL_*", intern=TRUE)


df_class = read.table(files[1], row.names=1, header=TRUE) %>%
				rownames_to_column("Gene") %>%
				as_tibble

res = apply(df_class[,-1], 2, function(x) which(x < 0.05))

png("~/www/test.png", height=1000, width=700)
upset(fromList(res), order.by = "freq")
dev.off()


# hist(apply(df_class, 1, max))






