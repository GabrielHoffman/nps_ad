

chr=$1

distance=$2

cell1=$3

cell2=$4

cell="${cell1}_${cell2}"

if [ ! -d "${cell1}_${cell2}" ]
then

mkdir ${cell1}_${cell2}

fi

#exit

export chr

export distance





gunzip -c /hpc/users/zengb02/Microglia_space/atacseq/zengb02/Jaro_caQTL_detection/meta-analysis/${cell2}/analysis/combined_results_chr${chr}_sorted.gz   >${cell}/merge_combined_results_for_chr${chr}_temp


awk '{print "chr"$1,$2,$3,$4}'  ${cell}/merge_combined_results_for_chr${chr}_temp  >${cell}/merge_combined_results_for_chr${chr}_temp_temp 



awk 'NR==FNR{b=$2" "$3;a[$4]=b;next} NR>FNR{print a[$4],$0}'   ${cell}/merge_combined_results_for_chr${chr}_temp_temp   ${cell}/merge_combined_results_for_chr${chr}_temp |perl -ne 'if(/^\S/){chomp;my @array=split;$array[3]=$array[0];$array[4]=$array[1];shift @array; shift @array;my $str=join("\t",@array);print "$str\n";}'  >${cell}/merge_combined_results_for_chr${chr}

rm ${cell}/*_chr${chr}_temp* 


#exit
module load bedtools 




cat  /hpc/users/zengb02/Microglia_space/atacseq/zengb02/Jaro_caQTL_detection/data/phenotype/AMPAD/Data_peaks_${cell2}.bed  |perl -ne '$_=~s/^chr//g;chomp;my @array=split;;my $start=$array[1]-$ENV{distance};if($start<0){$start=0;} my $end=$array[2]+$ENV{distance};print "$array[0]\t$start\t$end\t$array[-1]\n";' |sort -k1,1 -k2,2n |intersectBed   -a  ${cell}/merge_combined_results_for_chr${chr}   -b  - -wo |perl -ne 'chomp;my @array=split;if(!($array[6] eq $array[13])){next;} print "@array[0..9]\n";' |perl -ne 'chomp;my @array=split;my $str=join("\t",@array);print "$str\n";'>${cell}/adjust_merge_combined_results_for_chr${chr}_caQTL 

rm  ${cell}/merge_combined_results_for_chr${chr}










module load samtools

module load htslib


gunzip -c  /hpc/users/zengb02/CommonMind_arion/eQTL_catalogue/metabrain/2021-07-23/another/2021-07-23-cortex-EUR-80PCs-chr${chr}_converted.gz  |intersectBed   -a   ${cell}/adjust_merge_combined_results_for_chr${chr}_caQTL  -b  -  -wo  >${cell}/merge_eQTL_caQTL_chr${chr} 

cat ${cell}/merge_eQTL_caQTL_chr${chr} |perl -ne 'chomp;my @array=split;print "$array[6]\t$array[3]\t$array[8]\t$array[9]\t$array[16]\t$array[13]\t$array[18]\t$array[19]\n";' |awk '$2==$6'  |awk '$3!=0 && $7!=0' >${cell}/adjust_merge_eQTL_caQTL_chr${chr}

rm ${cell}/adjust_merge_combined_results_for_chr${chr}_caQTL  
rm ${cell}/merge_eQTL_caQTL_chr${chr} 


cut -f1 ${cell}/adjust_merge_eQTL_caQTL_chr${chr}   |sort |uniq  >${cell}/probe_list_chr${chr}


cut -f5 ${cell}/adjust_merge_eQTL_caQTL_chr${chr}   |sort |uniq  >${cell}/gene_list_chr${chr}


exit


cat   /hpc/users/zengb02/Microglia_space/atacseq/zengb02/Jaro_caQTL_detection/coloc_analysis/run_with_coloc/${cell}/adjust_merge_eQTL_caQTL_chr${chr}   |perl  -ne  'chomp;my @array=split;print "$array[0]\t$array[4]\n";' |sort  |uniq  >probe_gene_pair_${cell}_chr${chr}

split  -l 800  probe_gene_pair_${cell}_chr${chr}   ${cell}_chr${chr}_BIAO


export cell
export chr

perl -e  'my @file=glob("$ENV{cell}_chr$ENV{chr}_BIAO*");foreach my $file  (@file){system  " create_bsub_job.pl   --command    \" ./run_coloc.sh  $ENV{cell} $ENV{chr}  $file \" --time   5  --memory  11000  --output $ENV{cell}_$ENV{chr}_$file  --queue premium"; }'


#R --vanilla --slave --input /hpc/users/zengb02/Microglia_space/atacseq/zengb02/Jaro_caQTL_detection/coloc_analysis/Franke_eQTL/run_with_coloc/adjust_merge_eQTL_caQTL_chr${chr} --pair  < ./run_coloc.R


#perl  /hpc/users/zengb02/Microglia_space/atacseq/zengb02/Jaro_caQTL_detection/coloc_analysis/with_8_brain_cell_type_eQTL/another/estimate_coloc_v4_BZ.pl   ${cell}/adjust_merge_eQTL_caQTL_chr${chr}   >${cell}/adjust_adjust_merge_eQTL_caQTL_chr${chr}

#rm  ${cell}/adjust_merge_eQTL_caQTL_chr${chr}

