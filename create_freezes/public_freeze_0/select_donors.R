
# Dec 22, 2022
# Select cohort based on AD and controls

library(tidyverse)

# Load metdata
df_meta = read_csv("/sc/arion/projects/psychAD/NPS-AD/freeze1_rc/metadata/syn26527784_latest.csv")
rownames(df_meta) = df_meta$SubID

lst_clinical = readRDS("/sc/arion/projects/psychAD/NPS-AD/freeze1_rc/metadata/clinical_metadata_sampleSets_latest.RDS")

df_meta$Dx_AD = NA
i = df_meta$SubID %in% lst_clinical$controls_neuropathological
df_meta$Dx_AD[i] = "Control" 
i = df_meta$SubID %in% lst_clinical$AD
df_meta$Dx_AD[i] = "AD"
df_meta$Dx_AD[with(df_meta, PD | DLBD)] = NA

i = with(df_meta, Age > 60 & Brain_bank == "MSSM" & !is.na(Dx_AD))

# xtabs(~ Dx_AD + Brain_bank, df_meta[i,])
# xtabs(~ Dx_AD + PD, df_meta[i,])
# xtabs(~ Dx_AD + DLBD, df_meta[i,])


df = df_meta[i,]

xtabs(~ Dx_AD, df)
xtabs(~ Dx_AD + Sex, df)

ggplot(df, aes(Dx_AD, Age)) +
	geom_boxplot()

id_ctrl = df$SubID[df$Dx_AD=="Control"][1:150]
id_avail = df$SubID[df$Dx_AD=="AD"]
id_selected = c()

for( id in id_ctrl){

	query.Age = df$Age[df$SubID == id]
	query.Sex = df$Sex[df$SubID == id]

	df_target = df[df$SubID %in% id_avail,c("SubID", "Age", "Sex")]
	df_target = df_target[df_target$Sex == query.Sex,]

	id_new = df_target$SubID[which.min(abs(df_target$Age - query.Age))]

	id_selected = append(id_selected, id_new)
	id_avail = setdiff(id_avail, id_new)
}

df2 = df[df$SubID %in% c( id_selected, id_ctrl),]



xtabs(~ Dx_AD, df2)
xtabs(~ Dx_AD + Sex, df2)


ggplot(df2, aes(Dx_AD, Age)) + geom_boxplot()

df_request = data.frame(SubID1 = id_selected, SubID2 = id_ctrl)



file = "/sc/arion/projects/psychAD/NPS-AD/public_release_0/SubID_request.tsv"
write.table(df_request, file=file, row.names=FALSE, quote=FALSE, sep="\t")





