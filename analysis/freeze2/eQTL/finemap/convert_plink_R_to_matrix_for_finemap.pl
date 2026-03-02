#!/usr/bin/perl 
open(FHIN,"<$ARGV[0]");
my $wa=<FHIN>;
my %ld;
my %snp;
while(<FHIN>)
 {
   chomp;my @array=split;
#   print "array2 is $array[2]\n";
#   print "array5 is $array[5]\n";
   $snp{$array[2]}=$array[1];
   $snp{$array[5]}=$array[4];
   $ld{$array[2]}{$array[5]}=$array[6];
   $ld{$array[5]}{$array[2]}=$array[6];
 }
print "\t";
foreach my $snp (sort {$snp{$a} <=> $snp{$b}} keys %snp)
 {
   print "$snp\t";
 }
print "\n";


foreach my $snp1 (sort {$snp{$a} <=> $snp{$b}} keys %snp)
 {
   print "$snp1\t";
   foreach my $snp2 (sort {$snp{$a} <=> $snp{$b}} keys %snp)
     {
       if($snp1 eq $snp2)
        {
          print "1\t"; 
          next;
        }
       else 
        {
          if(exists($ld{$snp1}{$snp2}))
           {
             my $val=$ld{$snp1}{$snp2};
             print "$val\t";
           }
          elsif(exists($ld{$snp2}{$snp1}))
           {
             my $val=$ld{$snp2}{$snp1};
             print "$val\t";
           }
          else
           {
             print "There is some wrong, exit\n";
             exit(1);
           } 
        }
     }
    print "\n";
 }
 
   

