




cell=$1
step=$2
export cell
export step
###For test
#perl -e 'for(my $i=1;$i<=11;$i++){system "create_bsub_job.pl   --command \"  /hpc/users/zengb02/Microglia_space/atacseq/zengb02/Jaro_caQTL_detection/coloc_analysis/Franke_eQTL/run_with_coloc/another/extract.sh  $i 50000  $ENV{cell} neuron \"  --time 8 --memory 11000 --output extract_for_$i\_$ENV{cell}_neuron  --queue premium ";}'
###For test

perl -e 'for(my $i=12;$i<=22;$i++){system "create_bsub_job.pl   --command \"  /hpc/users/zengb02/CommonMind_arion/single_cell_eQTL_for_Gabriel_02212024/eQTL_detection_on_pearsonsubtype_level/dynamic_state_eQTL/Aging/run_dynamic_eQTL_with_fixed_theta/EN/another/colocalization_analysis/another/extract.sh  $i 50000  $ENV{cell} neuron  $ENV{step} \"  --time 8 --memory 11000 --output extract_for_$i\_$ENV{cell}_neuron  --queue premium ";}'



perl -e 'for(my $i=12;$i<=22;$i++){system "create_bsub_job.pl   --command \"  /hpc/users/zengb02/Microglia_space/atacseq/zengb02/Jaro_caQTL_detection/coloc_analysis/run_with_coloc/another/extract.sh  $i 50000  $ENV{cell} glia $ENV{step}\"  --time 8 --memory 11000 --output extract_for_$i\_$ENV{cell}_glia  --queue premium ";}'


###For test
#perl -e 'for(my $i=1;$i<=11;$i++){system "create_bsub_job.pl   --command \"  /hpc/users/zengb02/Microglia_space/atacseq/zengb02/Jaro_caQTL_detection/coloc_analysis/Franke_eQTL/run_with_coloc/another/extract.sh  $i 50000  $ENV{cell} glia \"  --time 8  --memory 11000 --output extract_for_$i\_$ENV{cell}_glia  --queue premium ";}'
###For test


