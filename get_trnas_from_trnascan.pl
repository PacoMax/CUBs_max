#!/usr/bin/perl

=head1 get_trnas_from_trnascan v.1.0

=head1 Created

<2018 - fgonzale>

=head1 The program gets a table with the number of tRNAs.

The user has to write
./get_trnas_from_trnascan.pl -i {input.table} > {output.table}


=head1 PARAMETERS

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
   if($_ =~ /\sTTT\s/) { $count[0]++ }
	elsif($_ =~ /\sTTC\s/) { $count[1]++ }
	elsif($_ =~ /\sTTA\s/) { $count[2]++ }
	elsif($_ =~ /\sTTG\s/) { $count[3]++ }
	elsif($_ =~ /\sTCT\s/) { $count[4]++ }
	elsif($_ =~ /\sTCC\s/) { $count[5]++ }
	elsif($_ =~ /\sTCA\s/) { $count[6]++ }
	elsif($_ =~ /\sTCG\s/) { $count[7]++ }
	elsif($_ =~ /\sTAT\s/) { $count[8]++ }
	elsif($_ =~ /\sTAC\s/) { $count[9]++ }
	elsif($_ =~ /\sTAA\s/) { $count[10]++ }
	elsif($_ =~ /\sTAG\s/) { $count[11]++ }
	elsif($_ =~ /\sTGT\s/) { $count[12]++ }
	elsif($_ =~ /\sTGC\s/) { $count[13]++ }
	elsif($_ =~ /\sTGA\s/) { $count[14]++ }
	elsif($_ =~ /\sTGG\s/) { $count[15]++ }
	elsif($_ =~ /\sCTT\s/) { $count[16]++ }
	elsif($_ =~ /\sCTC\s/) { $count[17]++ }
	elsif($_ =~ /\sCTA\s/) { $count[18]++ }
	elsif($_ =~ /\sCTG\s/) { $count[19]++ }
	elsif($_ =~ /\sCCT\s/) { $count[20]++ }
	elsif($_ =~ /\sCCC\s/) { $count[21]++ }
	elsif($_ =~ /\sCCA\s/) { $count[22]++ }
	elsif($_ =~ /\sCCG\s/) { $count[23]++ }
	elsif($_ =~ /\sCAT\s/) { $count[24]++ }
	elsif($_ =~ /\sCAC\s/) { $count[25]++ }
	elsif($_ =~ /\sCAA\s/) { $count[26]++ }
	elsif($_ =~ /\sCAG\s/) { $count[27]++ }
	elsif($_ =~ /\sCGT\s/) { $count[28]++ }
	elsif($_ =~ /\sCGC\s/) { $count[29]++ }
	elsif($_ =~ /\sCGA\s/) { $count[30]++ }
	elsif($_ =~ /\sCGG\s/) { $count[31]++ }
	elsif($_ =~ /\sATT\s/) { $count[32]++ }
	elsif($_ =~ /\sATC\s/) { $count[33]++ }
	elsif($_ =~ /\sATA\s/) { $count[34]++ }
	elsif($_ =~ /\sATG\s/) { $count[35]++ }
	elsif($_ =~ /\sACT\s/) { $count[36]++ }
	elsif($_ =~ /\sACC\s/) { $count[37]++ }
	elsif($_ =~ /\sACA\s/) { $count[38]++ }
	elsif($_ =~ /\sACG\s/) { $count[39]++ }
	elsif($_ =~ /\sAAT\s/) { $count[40]++ }
	elsif($_ =~ /\sAAC\s/) { $count[41]++ }
	elsif($_ =~ /\sAAA\s/) { $count[42]++ }
	elsif($_ =~ /\sAAG\s/) { $count[43]++ }
	elsif($_ =~ /\sAGT\s/) { $count[44]++ }
	elsif($_ =~ /\sAGC\s/) { $count[45]++ }
	elsif($_ =~ /\sAGA\s/) { $count[46]++ }
	elsif($_ =~ /\sAGG\s/) { $count[47]++ }
	elsif($_ =~ /\sGTT\s/) { $count[48]++ }
	elsif($_ =~ /\sGTC\s/) { $count[49]++ }
	elsif($_ =~ /\sGTA\s/) { $count[50]++ }
	elsif($_ =~ /\sGTG\s/) { $count[51]++ }
	elsif($_ =~ /\sGCT\s/) { $count[52]++ }
	elsif($_ =~ /\sGCC\s/) { $count[53]++ }
	elsif($_ =~ /\sGCA\s/) { $count[54]++ }
	elsif($_ =~ /\sGCG\s/) { $count[55]++ }
	elsif($_ =~ /\sGAT\s/) { $count[56]++ }
	elsif($_ =~ /\sGAC\s/) { $count[57]++ }
	elsif($_ =~ /\sGAA\s/) { $count[58]++ }
	elsif($_ =~ /\sGAG\s/) { $count[59]++ }
	elsif($_ =~ /\sGGT\s/) { $count[60]++ }
	elsif($_ =~ /\sGGC\s/) { $count[61]++ }
	elsif($_ =~ /\sGGA\s/) { $count[62]++ }
	elsif($_ =~ /\sGGG\s/) { $count[63]++ }
}
    for(my $j = 0; $j < @count; $j++) {
	print "$count[$j]\n";
    }
close MIFICH;

sub PrintHelp {
	system "pod2text -c $0";
	exit()
}
