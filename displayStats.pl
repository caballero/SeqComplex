#!/usr/bin/perl -w
use strict;
use FindBin;
use lib $FindBin::Bin;
use SeqComplex;
use Getopt::Long;
use Data::Dumper;

my %methodMetaData = (
  "at"  => { "title" => "Percent AT",
             "x_units" => "AT%" },
  "gc"  => { "title" => "Percent GC",
             "x_units" => "GC%" },
  "ats" => { "title" => "AT Skew",
             "x_units" => "(A-T)/(A+T)" },
  "gcs" => { "title" => "GC Skew",
             "x_units" => "(G-C)/(G+C)" },
  "ce" =>  { "title" => "Shannon Entropy",
             "x_units" => "bits" },
  "cm" =>  { "title" => "Markov Model Complexity",
             "x_units" => "entropy" },
  "cl" =>  { "title" => "Linguisitic Complexity",
             "x_units" => "word fraction" },
  "ct" =>  { "title" => "Trifonov Complexity",
             "x_units" => "word fraction" },
  "cz" =>  { "title" => "Compression Factor ( gzip )",
             "x_units" => "compression factor" },
  "cwf" =>  { "title" => "Wootton and Federhen Complexity",
             "x_units" => "complexity" },
  "cpg" => { "title" => "Percent CpGs",
             "x_units" => "CpG%" } );
           
  

=head1 NAME 

displayStats.pl [-graph] -in <STATSFILE>  > outputfile

=head1 DESCRIPTION

This script read in the output of gatherStats.pl and either 
output a data table or a set of Google Charts HTML.

=cut

my $help      = undef;
my $graph     = undef;
my $graphBins = 100;
my $seq       = "";
my $datafile  = "";

usage()
    if (
         !GetOptions(
                      'help|h'  => \$help,
                      'in=s'    => \$datafile,
                      'graph|g' => \$graph
         )
    );

usage() if ( defined $help );
usage() if ( $datafile eq "" );

# MAIN
my %hdrData = runAllMethodsInline( "ACGTACGT" );
my @methods = sort keys( %hdrData );
my $data    = serializeIN( $datafile );
my %rawData = %{$data};

if ( defined $graph )
{

  # Autobin each data column
  my %methodRanges = ();
  foreach my $method ( keys( %rawData ) )
  {
    $methodRanges{$method}->{'min'} = -1;
    $methodRanges{$method}->{'max'} = -1;
    foreach my $file ( keys( %{ $rawData{$method} } ) )
    {
      foreach my $value ( @{ $rawData{$method}->{$file} } )
      {
        if (    $methodRanges{$method}->{'min'} == -1
             || $value < $methodRanges{$method}->{'min'} )
        {
          $methodRanges{$method}->{'min'} = $value;
        }
        if (    $methodRanges{$method}->{'max'} == -1
             || $value > $methodRanges{$method}->{'max'} )
        {
          $methodRanges{$method}->{'max'} = $value;
        }
      }
    }
  }

  my %graphBinData = ();
  my %binSizes     = ();
  foreach my $method ( keys( %rawData ) )
  {
    my $rangeLen =
        $methodRanges{$method}->{'max'} - $methodRanges{$method}->{'min'};
    my $binSize = $rangeLen / $graphBins;
    if ( $binSize == 0 )
    {
      $binSize = $methodRanges{$method}->{'min'} / $graphBins;
      $methodRanges{$method}->{'min'} =
          $methodRanges{$method}->{'min'} - ( $graphBins * $binSize );
    }
    $binSizes{$method} = $binSize;
    foreach my $file ( keys( %{ $rawData{$method} } ) )
    {
      foreach my $value ( @{ $rawData{$method}->{$file} } )
      {
        $value -= $methodRanges{$method}->{'min'};
        my $bin = 0;
        if ( $binSize != 0 )
        {
          $bin = sprintf( "%0.0f", $value / $binSize );
        }
        $graphBinData{$method}->{$file}->[ $bin ]++;
      }
    }
  }

  print "<HTML>\n";
  print "<HEAD>\n";
  print "<SCRIPT type=\"text/javascript\" src=\"https://www.gstatic.com/charts/loader.js\"></SCRIPT>\n";
  print "  <SCRIPT type=\"text/javascript\">\n";

  print "    google.charts.load('current', {packages: ['corechart', 'line']});\n";
  print "    google.charts.setOnLoadCallback(drawChart);\n";
  print "    function drawChart() {\n";

  foreach my $method ( keys( %graphBinData ) )
  {
    my $binSize = $binSizes{$method};
    my $methodPrefix = $method;
    my $wordSize = 0;
    if ( $methodPrefix =~ /^([a-zA-Z]+)(\d*)$/ )
    {
      $methodPrefix = $1;
      $wordSize = $2 if ( $2 ne "" );
    }
 
    print "      var data_$method = new google.visualization.DataTable();\n";
    print "      data_$method.addColumn('number','" . $methodMetaData{$methodPrefix}->{'x_units'}  . "');\n";
    foreach my $file ( sort keys( %{ $graphBinData{$method} } ) )
    {
      print "      data_$method.addColumn('number',\'$file\');\n";
    }
    print "      data_$method.addRows([\n";
    my $bin = $methodRanges{$method}->{'min'};
    if ( $binSize == 0 )
    {
      my $str = "[$bin, ";
      if ( $methodPrefix =~ /at|gc|cpg/ )
      {
        $str = "[" . ( $bin * 100 ) . ", ";
      }
      foreach my $file ( sort keys( %{ $graphBinData{$method} } ) )
      {
        $str .= $graphBinData{$method}->{$file}->[ 0 ] . ",";
      }
      $str =~ s/,$//;
      $str .= " ] ] );\n";
      print "$str";
    } else
    {
      my $str = "";
      for ( my $i = 0 ; $i <= $graphBins ; $i++ )
      {
        if ( $methodPrefix =~ /at|gc|cpg/ )
        {
          $str .= "[" . ( $bin * 100 ) . ", ";
        }else
        {
          $str .= "[$bin, ";
        }
        foreach my $file ( sort keys( %{ $graphBinData{$method} } ) )
        {
          my $value = $graphBinData{$method}->{$file}->[ $i ];
          $value = 0 if ( !defined $value );
          $str .= "$value,";
        }
        $str =~ s/,$//;
        $str .= "],";
        $bin += $binSize;
      }
      $str =~ s/,$//;
      print "$str ] );\n";
    }
    print "var options_$method = {\n";

    print "               chart: {\n";
    print  "                 title: '" . $methodMetaData{$methodPrefix}->{'title'};
    if ( $wordSize > 0 ) 
    {
      print " ( word size = $wordSize )";
    }
    print "'\n";
    print "                      },\n";

    print "                 width: 900,\n";
    print "                 height: 500,\n";
    print "                 series: { 0: {axis: 'bins'} },\n";
    print "                 axes: { y: { bins: {label: '1k bins, 500bp overlap'} } }\n";
    print "};\n";

    print "var chart_$method = new google.charts.Line(document.getElementById('chart_$method'));\n";
    print "chart_$method.draw(data_$method, options_$method);\n";
  }
  print "}\n";
  print "</SCRIPT>\n";
  print "</HEAD>\n";
  print "<BODY>\n";
  foreach my $method ( sort keys( %graphBinData ) )
  {
    print "  <DIV ID=\"chart_$method\"></div>\n";
  }
  print "</BODY>\n";
  print "</HTML>\n";
}

=head1 SUBROUTINES

=cut

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
sub serializeIN
{
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
sub serializeOUT
{
  my $object      = shift;
  my $fileName    = shift;
  my $data_dumper = new Data::Dumper( [ $object ] );
  $data_dumper->Purity( 1 )->Terse( 1 )->Deepcopy( 1 );
  open OUT, ">$fileName";
  print OUT $data_dumper->Dump();
  close OUT;
}

=head1 AUTHOR

Robert Hubley
Juan Caballero
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

