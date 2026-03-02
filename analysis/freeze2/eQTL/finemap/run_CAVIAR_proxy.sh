


input=$1

gene=$2


chr=$3
num=`wc -l < ${input}`



if [[ ${num} == 0 ]]; then

exit

fi


Base=`basename ${input}`



export num

export Base

export gene
export input

#perl -e 'for(my $i=1;$i<=$ENV{num};$i++){for(my $j=1;$j<=$ENV{num};$j++){if($i==$j){print "\t1";} else{print "\t0";}} print "\n";}' |perl -ne 'chomp;my @array=split;my $str=join("\t",@array);print "$str\n";' >${gene}_${Base}_LD

./construct_LD.sh  ${input}  ${gene}  ${Base}  ${chr}


module load R/4.3.0

cat  ${input}_filter |perl -ne 'chomp;my @array=split;$array[1]=abs($array[1]);print "$array[0]\t$array[1]\n";' | sort -k2,2nr  |head -n 1  |perl -ne 'chomp;my @array=split;if($array[1]>=2){system "/hpc/packages/minerva-common/caviar/2017-04-06/CAVIAR-C++/CAVIAR     -l $ENV{gene}\_$ENV{Base}\_LD  -z $ENV{input}_filter  -c  1  -o $ENV{gene}\_$ENV{Base}\_caviar";}'


awk 'NR==FNR{a[$1]=$0;next} NR>FNR{print a[$1]}' ${gene}\_${Base}\_caviar_post   ${gene}\_${Base}\_caviar_set  >${gene}\_${Base}\_finemapping_result


rm ${gene}_${Base}_LD
rm ${gene}\_${Base}\_caviar* 

rm ${input}
