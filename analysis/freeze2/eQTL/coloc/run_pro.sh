










for cell2 in "EN" "Micro" "IN" "Astro" "Oligo" "OPC"
#for cell2 in "Micro"
do

for cell1 in "alzBellenguez/AD3" "alz2noapoe/AD" "asd/ASD"  "bip/BIP" "mdd2/MDD2" "pd/PD" "sz3/SCZ3"
do

export cell1
export cell2



perl -e 'for(my $i=1;$i<=22;$i++){my $cell=$ENV{cell1};my @array=split(/\//,$cell);system "create_bsub_job.pl   --command \" ./run_analysis.sh  $ENV{cell1} $array[-1]  $ENV{cell2}   $i \"   --time 2 --memory 25000 --output $array[-1]\_$ENV{cell2}_extract_for_$i  --queue premium";}'



done


done









