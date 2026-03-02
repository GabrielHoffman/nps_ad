


egene=$1   #a list of gene names 

gene_annotation=$2   #BED format gene annotation


chr=$2



if [ ! -d "chr${chr}" ]
then
mkdir chr${chr}
fi

cd chr${chr}












cat   ${gene_annotation} |awk 'NR==FNR{a[$4]=$0;next} NR>FNR{print a[$1]}'  -  ${egene} |perl -ne '$_=~s/^chr//;print $_;' |awk -v var="${chr}"  '$1==var' |sort |uniq |split -l 120 - 




export chr

queue=$1


export queue




 perl -e 'my @file=glob("x??");foreach my $file (@file){system "create_bsub_job.pl  --command \"sh   ./run_finemap1.sh  $file  $ENV{chr}\"  --memory 10000 --time 6  --output  MMeQTL_output_$file  --queue $ENV{queue} ";}'




