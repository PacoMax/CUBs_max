#!/usr/bin/perl

=head1 get_table_protein v.2.0

=head1 Created

<2019 - fgonzale>

=head1 The program gets a table with attributes from a single line multifasta. 
FILE_ID	Protein_ID

The user has to write
./get_table_protein.pl -f {fasta.file} > {table.out}

To conver multiple line to single line use:

perl -pe '/^>/ ? print "\n" : chomp' input.fasta | tail -n +2 > output.fasta

=head1 PARAMETERS

=over 4

=item B<-f> input fasta file

=item B<help|h> Display this help and exit

=back

=cut
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;


my %opts = ();
GetOptions (\%opts,'f=s',
		'h|help');

if(($opts{'h'}) || (scalar(keys(%opts)) == 0)){
   &PrintHelp();
}
unless($opts{f} && -f $opts{f}){
    &PrintHelp();
}


print "FILE_ID\tProtein_ID\n";
my $file_id = $opts{f};
my $protein_id;
open(MIFICH,$opts{f}) or die "Unable to open $opts{f} !\n ";
while(<MIFICH>) {
   chomp;
    if ($_=~/^>(.+)\s+./){
		$protein_id=$1;
		print "$file_id\t$protein_id\n"
	}
}
close(MIFICH);


sub PrintHelp {
	system "pod2text -c $0";
	exit()
}
