


file=$1

gene=$2

Chr=$3

Start=$4


End=$5



if [[ $file == *".gz" ]]; then

module load htslib  

tabix  ${file} ${Chr}:${Start}-${End} |awk -v var=${gene} '$7==var' |perl -ne 'chomp;my @array=split;if($array[-1]!=0){print "$array[3]\t$array[-1]\n";}'  >${gene}_adjust

fi




sh ./run_CAVIAR_proxy.sh    ${gene}_adjust  ${gene}  ${Chr}
