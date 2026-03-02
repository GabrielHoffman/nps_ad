



cell1=$1
cell11=$2
cell2=$3

chr=$4








gunzip -c  ${cell2}_combined_results_chr${chr}_sorted.gz   |awk  '$9!=0' | awk -v var=${chr} '$1==var' |awk  'NR==FNR{a[$4]=$0;next} NR>FNR{print  a[$4],$0}'  ./${cell1}_GWAS_summary_hg38_converted     -  |perl -ne  'if(/^\S/){print $_;}' |perl -ne  'chomp;my @array=split;print "$array[-6]\t$array[-9]\t$array[-4]\t$array[-3]\t$array[-2]\t$array[-1]\t$array[6]\t$array[3]\t$array[8]\t$array[9]\t$array[10]\t$array[11]\n";' >merged_infor_${cell11}_${cell2}_chr${chr}

cat merged_infor_${cell11}_${cell2}_chr${chr}   |perl -ne  'chomp;my @array=split;print "$array[0]\t$array[6]\n";' |sort |uniq  >gene_trait_pairs_${cell11}_${cell2}_chr${chr}

module load R

R --vanilla --slave --input ./merged_infor_${cell11}_${cell2}_chr${chr}     --pair  gene_trait_pairs_${cell11}_${cell2}_chr${chr}  < ./run_coloc.R








