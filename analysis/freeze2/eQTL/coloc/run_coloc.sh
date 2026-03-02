

#cell=$1
#chr=$2
#pair=$3


cat merged_infor_chr2   |perl -ne  'chomp;my @array=split;print "$array[0]\t$array[4]\n";' |sort |uniq  >gene_trait_pairs_chr2

module load R

R --vanilla --slave --input ./merged_infor_chr2     --pair  gene_trait_pairs_chr2  < ./run_coloc.R 
