#!/usr/bin/perl

=head1 get_trnas_from_tmp
=head1 Created
<2020 - fgonzale>
=head1 Uhe program gets a table with the tmps of tRNAs.
Uhe user has to write
./prepare_tmp.pl.pl -i {input.table} > {output.table}
=head1 PARAMEUERS
=over 4
=item B<-i> input table
=item B<help|h> Display this help and exit
=back
=cut
use strict;
use warnings;
use Getopt::Long;


my %opts = (); #declaramos que la entrada se va a llamar opts
GetOptions (\%opts,'i=s', #s de string 
		'l=s',
		'h|help');

if(($opts{'h'}) || (scalar(keys(%opts)) == 0)){
   &PrintHelp();
}
unless($opts{i} && -f $opts{i}){
    &PrintHelp();
}

open(MIFICH,$opts{i}) or die "Unable to open $opts{i} !\n ";


my @count= (0) x 64;

while(<MIFICH>) {
   chomp;
   if($_ =~ /UUU\t(.*)/) { $count[0]+=$1 }
	elsif($_ =~ /UUC\t(.*)/) { $count[1]+=$1 }
	elsif($_ =~ /UUA\t(.*)/) { $count[2]+=$1 }
	elsif($_ =~ /UUG\t(.*)/) { $count[3]+=$1 }
	elsif($_ =~ /UCU\t(.*)/) { $count[4]+=$1 }
	elsif($_ =~ /UCC\t(.*)/) { $count[5]+=$1 }
	elsif($_ =~ /UCA\t(.*)/) { $count[6]+=$1 }
	elsif($_ =~ /UCG\t(.*)/) { $count[7]+=$1 }
	elsif($_ =~ /UAU\t(.*)/) { $count[8]+=$1 }
	elsif($_ =~ /UAC\t(.*)/) { $count[9]+=$1 }
	elsif($_ =~ /UAA\t(.*)/) { $count[10]+=$1 }
	elsif($_ =~ /UAG\t(.*)/) { $count[11]+=$1 }
	elsif($_ =~ /UGU\t(.*)/) { $count[12]+=$1 }
	elsif($_ =~ /UGC\t(.*)/) { $count[13]+=$1 }
	elsif($_ =~ /UGA\t(.*)/) { $count[14]+=$1 }
	elsif($_ =~ /UGG\t(.*)/) { $count[15]+=$1 }
	elsif($_ =~ /CUU\t(.*)/) { $count[16]+=$1 }
	elsif($_ =~ /CUC\t(.*)/) { $count[17]+=$1 }
	elsif($_ =~ /CUA\t(.*)/) { $count[18]+=$1 }
	elsif($_ =~ /CUG\t(.*)/) { $count[19]+=$1 }
	elsif($_ =~ /CCU\t(.*)/) { $count[20]+=$1 }
	elsif($_ =~ /CCC\t(.*)/) { $count[21]+=$1 }
	elsif($_ =~ /CCA\t(.*)/) { $count[22]+=$1 }
	elsif($_ =~ /CCG\t(.*)/) { $count[23]+=$1 }
	elsif($_ =~ /CAU\t(.*)/) { $count[24]+=$1 }
	elsif($_ =~ /CAC\t(.*)/) { $count[25]+=$1 }
	elsif($_ =~ /CAA\t(.*)/) { $count[26]+=$1 }
	elsif($_ =~ /CAG\t(.*)/) { $count[27]+=$1 }
	elsif($_ =~ /CGU\t(.*)/) { $count[28]+=$1 }
	elsif($_ =~ /CGC\t(.*)/) { $count[29]+=$1 }
	elsif($_ =~ /CGA\t(.*)/) { $count[30]+=$1 }
	elsif($_ =~ /CGG\t(.*)/) { $count[31]+=$1 }
	elsif($_ =~ /AUU\t(.*)/) { $count[32]+=$1 }
	elsif($_ =~ /AUC\t(.*)/) { $count[33]+=$1 }
	elsif($_ =~ /AUA\t(.*)/) { $count[34]+=$1 }
	elsif($_ =~ /AUG\t(.*)/) { $count[35]+=$1 }
	elsif($_ =~ /ACU\t(.*)/) { $count[36]+=$1 }
	elsif($_ =~ /ACC\t(.*)/) { $count[37]+=$1 }
	elsif($_ =~ /ACA\t(.*)/) { $count[38]+=$1 }
	elsif($_ =~ /ACG\t(.*)/) { $count[39]+=$1 }
	elsif($_ =~ /AAU\t(.*)/) { $count[40]+=$1 }
	elsif($_ =~ /AAC\t(.*)/) { $count[41]+=$1 }
	elsif($_ =~ /AAA\t(.*)/) { $count[42]+=$1 }
	elsif($_ =~ /AAG\t(.*)/) { $count[43]+=$1 }
	elsif($_ =~ /AGU\t(.*)/) { $count[44]+=$1 }
	elsif($_ =~ /AGC\t(.*)/) { $count[45]+=$1 }
	elsif($_ =~ /AGA\t(.*)/) { $count[46]+=$1 }
	elsif($_ =~ /AGG\t(.*)/) { $count[47]+=$1 }
	elsif($_ =~ /GUU\t(.*)/) { $count[48]+=$1 }
	elsif($_ =~ /GUC\t(.*)/) { $count[49]+=$1 }
	elsif($_ =~ /GUA\t(.*)/) { $count[50]+=$1 }
	elsif($_ =~ /GUG\t(.*)/) { $count[51]+=$1 }
	elsif($_ =~ /GCU\t(.*)/) { $count[52]+=$1 }
	elsif($_ =~ /GCC\t(.*)/) { $count[53]+=$1 }
	elsif($_ =~ /GCA\t(.*)/) { $count[54]+=$1 }
	elsif($_ =~ /GCG\t(.*)/) { $count[55]+=$1 }
	elsif($_ =~ /GAU\t(.*)/) { $count[56]+=$1 }
	elsif($_ =~ /GAC\t(.*)/) { $count[57]+=$1 }
	elsif($_ =~ /GAA\t(.*)/) { $count[58]+=$1 }
	elsif($_ =~ /GAG\t(.*)/) { $count[59]+=$1 }
	elsif($_ =~ /GGU\t(.*)/) { $count[60]+=$1 }
	elsif($_ =~ /GGC\t(.*)/) { $count[61]+=$1 }
	elsif($_ =~ /GGA\t(.*)/) { $count[62]+=$1 }
	elsif($_ =~ /GGG\t(.*)/) { $count[63]+=$1 }
}
    for(my $j = 0; $j < @count; $j++) {
	print "$count[$j]\n";
    }
close MIFICH;

sub PrintHelp {
	system "pod2text -c $0";
	exit()
}
