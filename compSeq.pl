#!/usr/bin/perl -w
use strict;
use SeqComplex;
use Getopt::Long;

=head1 NAME 

compSeq.pl -s STEP -win WINDOW -f FASTA

=head1 DESCRIPTION

This script read a big fasta sequence and calculate a collection of 
composition and complexity methods. We return a big data matrix. 

=cut

my $help    = undef;
my $file    = undef;
my $step    = 500;
my $win     = 1000;
my $nmin    = 0.3;
my $seq     = '';
my @methods = qw/gc at gcs ats cpg cwf ce cz cm1 cm2 cm3 cm4 cm5 cm6 ct1 ct2 ct3 ct4 ct5 ct6 cl1 cl2 cl3 cl4 cl5 cl6/;


usage() if (!GetOptions( 
				'help|h'     => \$help,
				'fasta|f=s'  => \$file,
				'step|s:i'   => \$step,
				'win|w:i'    => \$win,
				'nmin|n:i'   => \$nmin
			));
usage() if (defined $help);
usage() unless (defined $file);	

my $size = readFasta($file);
my $end  = $size - $win; 

# MAIN
#print "#POSITION\t", join("\t", @methods), "\n";
my $chr = $file;
$chr =~ s/\.fa//;

for (my $pos = 0; $pos <= $end; $pos += $step) {
	my $str     = substr($seq, $pos, $win);
	my $nnum    = $str =~ tr/N/N/;
	next if ($nnum / $win > $nmin);
	my %results = runAllMethods($str);
	my $a       = $pos + 1;
	my $b       = $a + $win;
	print join ("\t", $chr, $a, $b);
	foreach my $m (@methods) {
		my $val = shift @{ $results{$m} };
		print "\t$val";
	}
	print "\n";
}

=head1 SUBROUTINES

=cut

sub readFasta {
	my $f = shift @_;
	open F, "$f" or die "cannot open $f\n";
	while (<F>) {
		next if(/>/);
		chomp;
		$seq .= uc($_);
	}
	close F;
	return length $seq;
}

sub usage {
print <<_HERE_
Usage: compSeq.pl -s STEP -win WINDOW -f FASTA

DEFAULTS:
  step    = $step   Overlap window
  window  = $win    Fragment size
  nmin    = $nmin   Minimal N's fraction
_HERE_
;
	exit 1;
}

=head1 AUTHOR

Juan Caballero
Institute for Systems Biology @ 2010

=head1 CONTACT

jcaballero@systemsbiology.org

=head1 LICENSE

This is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with code.  If not, see <http://www.gnu.org/licenses/>.

=cut
