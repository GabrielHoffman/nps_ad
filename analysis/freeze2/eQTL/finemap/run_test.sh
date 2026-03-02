

input=$1


chr=$2

data=~/CommonMind_arion/single_cell_eQTL_for_Gabriel_02212024/data/genotype/genotype_for_sceQTL_04032023_chr${chr}_geno_filter_maf_filter_with_both_genotype_and_phenotype

cut -f1  ${input}  |sort |uniq  >variants.txt

#plink --bfile ${data} --extract variants.txt --r2  --out ld


plink --bfile ${data} --extract variants.txt --r  --ld-window-r2   0   --ld-window-kb 5000  --ld-window 99999   --out ld

#R CMD BATCH  ./contruct_LD.R     ld.ld   ${gene}_${Base}_LD 


perl   convert_plink_R_to_matrix_for_finemap.pl    ld.ld >USP37_test_LD 




#Rscript  ./contruct_LD.R  ${input}  ld.ld  ${gene}_${Base}_LD
