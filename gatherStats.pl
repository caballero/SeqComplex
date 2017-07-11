#!/usr/bin/perl -w
use strict;
use FindBin;
use lib $FindBin::Bin;
use SeqComplex;
use Getopt::Long;
use Data::Dumper;

=head1 NAME 

gatherStats.pl -s STEP -win WINDOW -out <OUTFILE> FASTA <FASTA ..> 

=head1 DESCRIPTION

This script read a big fasta sequence and calculate a collection of 
composition and complexity methods. We return a big data matrix. 

=cut

my $help      = undef;
my $baseline  = undef;
my $step      = 500;
my $win       = 1000;
my $nmin      = 0.3;
my $graph     = undef;
my $graphBins = 100;
my $seq       = "";
my $outfile   = "";

usage()
    if (
         !GetOptions(
                      'help|h'     => \$help,
                      'out:s'      => \$outfile,
                      'step|s:i'   => \$step,
                      'win|w:i'    => \$win,
                      'nmin|n:i'   => \$nmin
         )
    );

usage() if ( defined $help );
usage() unless ( defined @ARGV && $outfile ne "" );

# MAIN
my %hdrData = runAllMethodsInline( "ACGTACGT" );
my @methods = sort keys( %hdrData );
my %rawData = ();

foreach my $file ( @ARGV )
{
  $seq  = '';
  my $size = readFasta( $file );
  my $end  = $size - $win;
  for ( my $pos = 0 ; $pos <= $end ; $pos += $step )
  {
    my $str = substr( $seq, $pos, $win );
    my $nnum = $str =~ tr/N/N/;
    next if ( $nnum / $win > $nmin );
    my %results = runAllMethodsInline( $str );
    foreach my $m ( @methods )
    {
      my $val = shift @{ $results{$m} };
      push @{ $rawData{$m}->{$file} }, $val;
    }
  }
}
serializeOUT( \%rawData, $outfile );


=head1 SUBROUTINES

=cut

sub readFasta
{
  my $f = shift @_;
  open F, "$f" or die "cannot open $f\n";
  while ( <F> )
  {
    next if ( />/ );
    chomp;
    $seq .= uc( $_ );
  }
  close F;
  return length $seq;
}

sub usage
{
  exec "pod2text $0";
  exit 1;
}

##-------------------------------------------------------------------------##
## Use: my $val = serializeIN( $filename );
##
##      $filename       : A filename containing a serialized object
##
##  Returns
##
##      Uses the Data::Dumper module to read in data
##      from a serialized PERL object or data structure.
##
##-------------------------------------------------------------------------##
sub serializeIN {
  my $fileName     = shift;
  my $fileContents = "";
  my $oldSep       = $/;
  undef $/;
  my $in;
  open $in, "$fileName";
  $fileContents = <$in>;
  $/            = $oldSep;
  close $in;
  return eval( $fileContents );
}

##-------------------------------------------------------------------------##
## Use: serializeOUT( $data, $filename );
##
##        $filename     : A filename to be created
##
##  Returns
##
##      Uses the Data::Dumper module to save out the data
##      structure as a text file.  This text file can be
##      read back into an object of this type.
##
##-------------------------------------------------------------------------##
sub serializeOUT {
  my $object      = shift;
  my $fileName    = shift;
  my $data_dumper = new Data::Dumper( [ $object ] );
  $data_dumper->Purity( 1 )->Terse( 1 )->Deepcopy( 1 );
  open OUT, ">$fileName";
  print OUT $data_dumper->Dump();
  close OUT;
}


=head1 AUTHOR

Juan Caballero
Robert Hubley
Institute for Systems Biology @ 2010-2015

=head1 CONTACT

rhubley@systemsbiology.org

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

