



file=$1

chr=$2

export chr

cat ${file} |perl -ne 'chomp;my @array=split;$array[0]=~s/^chr//g;my $start=$array[1]-1000000;my $end=$array[2]+1000000;if($start<0){$start=0;} system " ./finemap1.sh  eQTL_summary_results_chr$ENV{chr}_sorted.gz  $array[-1]  $array[0] $start $end";'

