





cat chr*/*_finemapping_result  |perl -ne 'chomp;my @array=split;print "$array[0]\t$array[-1]\n";' >merged_finemapping_results
