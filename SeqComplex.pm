package SeqComplex;

=head1 NAME

SeqComplex

=head1 SYNOPSIS

  # Example 1. Run all methods for each sequence
  use SeqComplex;

  my @seq = loadSeqs($file);
  my $num = 0;
  foreach my $seq (@seq) {
	$num++;
	my %results = Complex::runAllMethods( $seq, $win, $kmer );
	foreach my $m (keys %results) {
		my $res = join "\t", @{ $results{$m} };
		print "seq_$num\t$m\t$res\n";
	}
  }

  ## Example 2. Run a particular method
  ## NOTE: all methods require references as inputs and returns an ARRAY.
  use SeqComplex qw(:methods)

  # Calculate the GC content every 1kb
  my $seq = 'ATGC' x 10000;
  my $win = 1000;
  my @gc_vals = gc(\$seq, \$win);

=head1 DESCRIPTION

Calculate composition and complexity of a DNA sequence.

=head2 EXPORT

runAllMethods

=head2 EXPORT_OK

:methods exports methods  = [gc at gcs ats cpg cwf ce cz ct cl cm]
:utils   exports basic utils = [log_k pot createWords countWords randSeq]
:all     exports @methods and @utils

=cut

use 5.008008;
use strict;
use warnings;
use Compress::Zlib;
use Data::Dumper;

require Exporter;

our @ISA = qw(Exporter);

our @methods = qw/gc at gcs ats cpg cwf ce cz ct cl cm/;

our @utils = qw/log_k pot createWords countWords randSeq/;

our %EXPORT_TAGS = (
                     'all'     => [ @methods, @utils ],
                     'methods' => [ @methods ],
                     'utils'   => [ @utils ]
);

our @EXPORT_OK = (
                   @{ $EXPORT_TAGS{'all'} },
                   @{ $EXPORT_TAGS{'methods'} },
                   @{ $EXPORT_TAGS{'utils'} }
);

our @EXPORT = qw( runAllMethods runAllMethodsInline );

our $VERSION = '0.01';

our $k = 4;    # alphabet size <A, C, G, T>

our @alphabet = qw/A G C T/;

=head1 SUBROUTINES

=head2 Methods

=cut

=head3 runAllMethods

Function to run each method in a DNA sequence.

Call: runAllMethods ( $seq, $win, $word ) STRING, NUMBER, NUMBER

Return: %values HASH (KEY = method, VALUES = @values)

=cut

sub runAllMethods
{
  my $gseq = shift @_;
  my $gwin = shift @_;
  $gwin ||= length $gseq;    # Use full length of the sequence
  my $gword  = shift @_;
  my @gwords = ();
  if ( defined $gword ) { push @gwords, $gword; }
  else { @gwords = ( 1, 2, 3, 4, 5, 6 ); }

  my %gvalues = ();
  @{ $gvalues{'gc'} }  = gc( \$gseq,  \$gwin );
  @{ $gvalues{'gcs'} } = gcs( \$gseq, \$gwin );
  @{ $gvalues{'cpg'} } = cpg( \$gseq, \$gwin );
  @{ $gvalues{'ce'} }  = ce( \$gseq,  \$gwin );
  @{ $gvalues{'cz'} }  = cz( \$gseq,  \$gwin );
  @{ $gvalues{'cwf'} } = cwf( \$gseq, \$gwin );

  #	@{ $gvalues{'at' } } =  at( \$gseq, \$gwin );
  @{ $gvalues{'ats'} } = ats( \$gseq, \$gwin );

  #	@{ $gvalues{'ket'} } = ket( \$gseq, \$gwin );
  #	@{ $gvalues{'pur'} } = pur( \$gseq, \$gwin );
  foreach my $ws ( @gwords )
  {
    @{ $gvalues{"ct$ws"} } = ct( \$gseq, \$gwin, \$ws );
    @{ $gvalues{"cl$ws"} } = cl( \$gseq, \$gwin, \$ws );
    @{ $gvalues{"cm$ws"} } = cm( \$gseq, \$gwin, \$ws );
  }
  return %gvalues;
}

=head3 runAllMethodsInline

Function to run each method in a DNA sequence -- all methods inlined

Call: runAllMethodsInline ( $seq, $win, $word ) STRING, NUMBER, NUMBER

Return: %values HASH (KEY = method, VALUES = @values)

=cut
sub runAllMethodsInline
{
  my $seq = shift @_;
  my $win = shift @_;
  $win ||= length $seq;    # Use full length of the sequence
  my $gword  = shift @_;
  my @gwords = ();

  if ( defined $gword ) { push @gwords, $gword; }
  else { @gwords = ( 1, 2, 3, 4, 5, 6 ); }
  my @sortedGWords = sort { $b <=> $a } @gwords;
  my $len = length $seq;
  my %gvalues;
  my $tot  = 0;
  my $cnt1 = 0;
  my $cnt2 = 0;
  my $val  = 0;
  my $up   = 0;

  $up = log_k( $k, $win ) * $win;
  for ( my $p = 0 ; $p <= ( $len - $win ) ; $p += $win )
  {
    my %pots = ();
    my %words = ();

    # Init state when word == 1
    foreach my $b ( @alphabet )
    {
      $words{1}{$b} = 0;
    }

    my $str   = substr( $seq, $p, $win );
    my $containsNonBase = 0;
    for ( my $i = 0 ; $i < $win ; $i++ )
    {
      my $rem = $win - $i > $sortedGWords[ 0 ] ? $sortedGWords[ 0 ] : $win - $i;
      my $tmpSeq = substr( $str, $i, $rem );
      if ( $tmpSeq !~ /[^ACGT]/ )
      {
        $words{$rem}{$tmpSeq}++;
        $containsNonBase = 0;
      }else
      {
        $containsNonBase = 1;
      }
      for ( my $j = 1 ; $j <= $#sortedGWords ; $j++ )
      {
        next if ( $sortedGWords[ $j ] >= $rem );
        my $elem = substr( $tmpSeq, 0, $sortedGWords[ $j ] );
        next if ( $containsNonBase && $elem =~ /[^ACGT]/ );
        $words{ $sortedGWords[ $j ] }{$elem}++;
      }
    }

    foreach my $ws ( @gwords )
    {
      $pots{$ws} = pot( $k, $ws );
    }
    my %elm    = %{ $words{1} };
    my $totNuc = $elm{'C'} + $elm{'G'} + $elm{'T'} + $elm{'A'};

    # GC
    $val = ( $elm{'C'} + $elm{'G'} ) / $totNuc if ( $totNuc > 1 );
    push @{ $gvalues{'gc'} }, $val;

    # AT
    $val = ( $elm{'A'} + $elm{'T'} ) / $totNuc if ( $totNuc > 1 );
    push @{ $gvalues{'at'} }, $val;

    # GCS
    $cnt1 = $elm{'G'} + $elm{'C'};
    $val = abs( $elm{'G'} - $elm{'C'} ) / $cnt1 if ( $cnt1 > 1 );
    push @{ $gvalues{'gcs'} }, $val;

    # ATS
    $cnt1 = $elm{'A'} + $elm{'T'};
    $val = abs( $elm{'A'} - $elm{'T'} ) / $cnt1 if ( $cnt1 > 1 );
    push @{ $gvalues{'ats'} }, $val;

    # CPG
    $cnt1 = $str =~ tr/CG/CG/;
    $cnt2 = $elm{'G'} * $elm{'C'};
    $val  = $cnt1 / $cnt2 if ( $cnt2 > 1 );
    push @{ $gvalues{'cpg'} }, $val;

    # CWF
    my $dw = 0;
    foreach my $b ( keys %elm )
    {
      next unless ( $elm{$b} > 0 );
      $dw += log_k( $k, $elm{$b} );
    }
    $val = ( $up - $dw ) / $totNuc if ( $totNuc > 1 );
    push @{ $gvalues{'cwf'} }, $val;

    # CE
    $val = 0;
    foreach my $b ( keys %elm )
    {
      next unless ( $elm{$b} > 0 );
      my $r = $elm{$b} / $totNuc if ( $totNuc > 1 );
      $val -= $r * log_k( $k, $r );
    }
    push @{ $gvalues{"ce"} }, $val;

    # CM
    foreach my $ws ( @gwords )
    {
      my $cm = 0;
      my $dw = $win - $ws - 1;
      foreach my $b ( keys %{ $words{$ws} } )
      {
        next unless ( $words{$ws}{$b} > 0 );
        my $r = $words{$ws}{$b} / $dw;
        $cm -= $r * log_k( $pots{$ws}, $r );
      }
      push @{ $gvalues{"cm$ws"} }, $cm;
    }

    # CL
    foreach my $ws ( @gwords )
    {
      my $sum_vl = 0;
      my $sum_vm = 0;
      for ( my $l = 1 ; $l <= $ws ; $l++ )
      {
        my $vl  = 0;
        my $vm  = 0;
        my $pot = $pots{$l};
        if ( $pot < ( $win - $l + 1 ) ) { $vm = $pot; }
        else { $vm = $win - $l + 1; }
        $sum_vm += $vm;
        foreach my $b ( keys %{ $words{$l} } )
        {
          next unless ( $words{$l}{$b} > 0 );
          $vl++;
        }
        $sum_vl += $vl;
      }
      my $r = 0;
      $r = $sum_vl / $sum_vm if ( $sum_vm > 0 );
      push @{ $gvalues{"cl$ws"} }, $r;
    }

    # CT
    foreach my $ws ( @gwords )
    {
      my $ct = 1;
      for ( my $l = 1 ; $l <= $ws ; $l++ )
      {
        my $vl  = 0;
        my $vm  = 0;
        my $pot = $pots{$l};
        if ( $pot < ( $win - $l + 1 ) ) { $vm = $pot; }
        else { $vm = $win - $l + 1; }
        foreach my $b ( keys %{ $words{$l} } )
        {
          next unless ( $words{$l}{$b} > 0 );
          $vl++;
        }
        $ct *= $vl / $vm;
      }
      push @{ $gvalues{"ct$ws"} }, $ct;
    }

    # CZ
    my $r   = 0;
    my $tmp = 'temp' . int( rand 1e8 );
    open my $th, ">", $tmp or die "cannot create temp file\n";
    print $th $str;
    close $th;
    my $size = -s $tmp;
    system( "gzip $tmp" );
    $tmp .= ".gz";

    if ( -e $tmp and -s $tmp )
    {
      my $cmpz = -s $tmp;
      $r = $size / $cmpz;
    } else
    {
      $r = 'NA';
    }
    unlink $tmp;
    push @{ $gvalues{"cz"} }, $r;
  }
  return %gvalues;
}

=head3 gc

Function to calculate the GC content.

Call: gc( \$seq, \$win ) STRING, NUMBER

Return: @values ARRAY

=cut

sub gc
{
  my $seq    = shift;
  my $win    = shift;
  my $len    = length $$seq;
  my @values = ();
  for ( my $p = 0 ; $p <= ( $len - $$win ) ; $p += $$win )
  {
    my $str = substr( $$seq, $p, $$win );
    my %elm = countWords( $str, 1 );
    my $r   = 0;
    my $tot = $elm{'C'} + $elm{'G'} + $elm{'T'} + $elm{'A'};
    $r = ( $elm{'C'} + $elm{'G'} ) / $tot if ( $tot > 1 );
    push @values, $r;
  }
  return @values;
}

=head3 at

Function to calculate the AT content.

Call: at( \$seq, \$win ) STRING, NUMBER

Return: @values ARRAY

=cut

sub at
{
  my $seq    = shift;
  my $win    = shift;
  my $len    = length $$seq;
  my @values = ();
  for ( my $p = 0 ; $p <= ( $len - $$win ) ; $p += $$win )
  {
    my $str = substr( $$seq, $p, $$win );
    my %elm = countWords( $str, 1 );
    my $r   = 0;
    my $tot = $elm{'C'} + $elm{'G'} + $elm{'T'} + $elm{'A'};
    $r = ( $elm{'A'} + $elm{'T'} ) / $tot if ( $tot > 1 );
    push @values, $r;
  }
  return @values;
}

=head3 gcs

Function to calculate the GC skew content.

Call: gcs( \$seq, \$win ) STRING, NUMBER

Return: @values ARRAY

=cut

sub gcs
{
  my $seq    = shift;
  my $win    = shift;
  my $len    = length $$seq;
  my @values = ();
  for ( my $p = 0 ; $p <= ( $len - $$win ) ; $p += $$win )
  {
    my $str = substr( $$seq, $p, $$win );
    my %elm = countWords( $str, 1 );
    my $l   = $elm{'G'} + $elm{'C'};
    my $r   = 0;
    $r = abs( $elm{'G'} - $elm{'C'} ) / $l if ( $l > 1 );
    push @values, $r;
  }
  return @values;
}

=head3 cpg

Function to calculate the CpG content.

Call: cpg( \$seq, \$win ) STRING, NUMBER

Return: @values ARRAY

=cut

sub cpg
{
  my $seq    = shift;
  my $win    = shift;
  my $len    = length $$seq;
  my @values = ();
  for ( my $p = 0 ; $p <= ( $len - $$win ) ; $p += $$win )
  {
    my $str = substr( $$seq, $p, $$win );
    my %elm = countWords( $str, 1 );
    my $c = $str =~ tr/CG/CG/;
    my $l = $elm{'G'} * $elm{'C'};
    my $r = 0;
    $r = $c / $l if ( $l > 1 );
    push @values, $r;
  }
  return @values;
}

=head3 ats

Function to calculate the AT skew content.

Call: ats( \$seq, \$win ) STRING, NUMBER

Return: @values ARRAY

=cut

sub ats
{
  my $seq    = shift;
  my $win    = shift;
  my $len    = length $$seq;
  my @values = ();
  for ( my $p = 0 ; $p <= ( $len - $$win ) ; $p += $$win )
  {
    my $str = substr( $$seq, $p, $$win );
    my %elm = countWords( $str, 1 );
    my $l   = $elm{'A'} + $elm{'T'};
    my $r   = 0;
    $r = abs( $elm{'A'} - $elm{'T'} ) / $l if ( $l > 1 );
    push @values, $r;
  }
  return @values;
}

=head3 ket

Function to calculate the Keto skew content.

  G+C-A-T / Total

Call: ket( \$seq, \$win ) STRING, NUMBER

Return: @values ARRAY

=cut

sub ket
{
  my $seq    = shift;
  my $win    = shift;
  my $len    = length $$seq;
  my @values = ();
  for ( my $p = 0 ; $p <= ( $len - $$win ) ; $p += $$win )
  {
    my $str = substr( $$seq, $p, $$win );
    my %elm = countWords( $str, 1 );
    my $r   = 0;
    my $tot = $elm{'G'} + $elm{'T'} + $elm{'A'} + $elm{'C'};
    $r = abs( $elm{'G'} + $elm{'T'} - $elm{'A'} - $elm{'C'} ) / $tot
        if ( $tot > 1 );
    push @values, $r;
  }
  return @values;
}

=head3 pur

Function to calculate the Purine skew content.

Call: pur( \$seq, \$win ) STRING, NUMBER

Return: @values ARRAY

=cut

sub pur
{
  my $seq    = shift;
  my $win    = shift;
  my $len    = length $$seq;
  my @values = ();
  for ( my $p = 0 ; $p <= ( $len - $$win ) ; $p += $$win )
  {
    my $str = substr( $$seq, $p, $$win );
    my %elm = countWords( $str, 1 );
    my $r   = 0;
    my $tot = $elm{'G'} + $elm{'T'} + $elm{'A'} + $elm{'C'};
    $r = abs( $elm{'G'} - $elm{'T'} + $elm{'A'} - $elm{'C'} ) / $tot
        if ( $tot > 1 );
    push @values, $r;
  }
  return @values;
}

=head3 cwf

Function to calculate the Complexity by Wootton & Federhen values.

Call: cwf( \$seq, \$win, ) STRING, NUMBER

Return: @values ARRAY

=cut

sub cwf
{
  my $seq    = shift;
  my $win    = shift;
  my $len    = length $$seq;
  my @values = ();
  my $up     = 0;
  for ( my $i = 1 ; $i <= $$win ; $i++ ) { $up += log_k( $k, $$win ); }
  for ( my $p = 0 ; $p <= ( $len - $$win ) ; $p += $$win )
  {
    my $str = substr( $$seq, $p, $$win );
    my %elm = countWords( $str, 1 );
    my $r   = 0;
    my $dw  = 0;
    my $tot = $elm{'G'} + $elm{'T'} + $elm{'A'} + $elm{'C'};

    foreach my $b ( keys %elm )
    {
      next unless ( $elm{$b} > 0 );
      $dw += log_k( $k, $elm{$b} );
    }
    $r = ( $up - $dw ) / $tot if ( $tot > 1 );
    push @values, $r;
  }
  return @values;
}

=head3 ce

Function to calculate the Complexity Entropy values.

Call: ce( \$seq, \$win ) STRING, NUMBER

Return: @values ARRAY

=cut

sub ce
{
  my $seq    = shift;
  my $win    = shift;
  my $len    = length $$seq;
  my @values = ();
  for ( my $p = 0 ; $p <= ( $len - $$win ) ; $p += $$win )
  {
    my $str = substr( $$seq, $p, $$win );
    my %elm = countWords( $str, 1 );
    my $ce  = 0;
    my $tot = $elm{'G'} + $elm{'T'} + $elm{'A'} + $elm{'C'};

    foreach my $b ( keys %elm )
    {
      next unless ( $elm{$b} > 0 );
      my $r = 0;
      $r = $elm{$b} / $tot if ( $tot > 1 );
      $ce -= $r * log_k( $k, $r );
    }
    push @values, $ce;
  }
  return @values;
}

=head3 cm

Function to calculate the Complexity in Markov model values.

Call: cm( \$seq, \$win, \$word ) STRING, NUMBER, NUMBER

Return: @values ARRAY

=cut

sub cm
{
  my $seq    = shift;
  my $win    = shift;
  my $word   = shift;
  my $len    = length $$seq;
  my @values = ();
  my $m      = pot( $k, $$word );
  for ( my $p = 0 ; $p <= ( $len - $$win ) ; $p += $$win )
  {
    my $str = substr( $$seq, $p, $$win );
    my %elm = countWords( $str, $$word );
    my $cm  = 0;
    my $dw  = $$win - $$word - 1;
    foreach my $b ( keys %elm )
    {
      next unless ( $elm{$b} > 0 );
      my $r = $elm{$b} / $dw;
      $cm -= $r * log_k( $m, $r );
    }
    push @values, $cm;
  }
  return @values;
}

=head3 cl

Function to calculate the Complexity Linguistic values.

Call: cl( \$seq, \$win, \$word ) STRING, NUMBER, NUMBER

Return: @values ARRAY

=cut

sub cl
{
  my $seq    = shift;
  my $win    = shift;
  my $word   = shift;
  my $len    = length $$seq;
  my @values = ();
  for ( my $p = 0 ; $p <= ( $len - $$win ) ; $p += $$win )
  {
    my $str    = substr( $$seq, $p, $$win );
    my $sum_vl = 0;
    my $sum_vm = 0;
    for ( my $l = 1 ; $l <= $$word ; $l++ )
    {
      my $vl  = 0;
      my $vm  = 0;
      my $pot = pot( $k, $l );
      if ( $pot < ( $$win - $l + 1 ) ) { $vm = $pot; }
      else { $vm = $$win - $l + 1; }
      $sum_vm += $vm;
      my %elm = countWords( $str, $l );
      foreach my $b ( keys %elm )
      {
        next unless ( $elm{$b} > 0 );
        $vl++;
      }
      $sum_vl += $vl;
    }
    my $r = 0;
    $r = $sum_vl / $sum_vm if ( $sum_vm > 0 );
    push @values, $r;
  }
  return @values;
}

=head3 ct

Function to calculate the Complexity by Trifonov values.

Call: ct( \$seq, \$win, \$word ) STRING, NUMBER, NUMBER

Return: @values ARRAY

=cut

sub ct
{
  my $seq    = shift;
  my $win    = shift;
  my $word   = shift;
  my $len    = length $$seq;
  my @values = ();
  for ( my $p = 0 ; $p <= ( $len - $$win ) ; $p += $$win )
  {
    my $ct = 1;
    my $str = substr( $$seq, $p, $$win );
    for ( my $l = 1 ; $l <= $$word ; $l++ )
    {
      my $vl  = 0;
      my $vm  = 0;
      my $pot = pot( $k, $l );
      if ( $pot < ( $$win - $l + 1 ) ) { $vm = $pot; }
      else { $vm = $$win - $l + 1; }
      my %elm = countWords( $str, $l );
      foreach my $b ( keys %elm )
      {
        next unless ( $elm{$b} > 0 );
        $vl++;
      }
      $ct *= $vl / $vm;
    }
    push @values, $ct;
  }
  return @values;
}

=head3 cz

Function to calculate the Compression factor.

Call: cz( \$seq, \$win ) STRING, NUMBER

Return: @values ARRAY

=cut

sub cz
{
  my $seq    = shift;
  my $win    = shift;
  my $len    = length $$seq;
  my @values = ();
  for ( my $p = 0 ; $p <= ( $len - $$win ) ; $p += $$win )
  {
    my $str = substr( $$seq, $p, $$win );
    my $r   = 0;
    my $tmp = 'temp' . int( rand 1e8 );
    open my $th, ">", $tmp or die "cannot create temp file\n";
    print $th $str;
    close $th;
    my $size = -s $tmp;
    system( "gzip $tmp" );
    $tmp .= ".gz";

    if ( -e $tmp and -s $tmp )
    {
      my $cmpz = -s $tmp;
      $r = $size / $cmpz;
    } else
    {
      $r = 'NA';
    }
    unlink $tmp;
    push @values, $r;
  }
  return @values;
}

=head3 clz

Function to calculate the CLZ values.

Call: clz( \$seq, \$win ) STRING, NUMBER

Return: @values ARRAY

* Not implemented yet!

=cut

sub clz
{
  print "Sorry, not implemented yet!\n";
  exit;
}

=head2 UTILITIES

=head3 pot

Function for calculate the exponential of a number.

Call: pot( $num, $exp ) NUMBER, NUMBER

Return: $res NUMBER

=cut

sub pot
{
  my $num = shift @_;
  my $exp = shift @_;
  if ( $num == 4 )
  {
    if    ( $exp == 1 )  { return 4; }
    elsif ( $exp == 2 )  { return 16; }
    elsif ( $exp == 3 )  { return 64; }
    elsif ( $exp == 4 )  { return 256; }
    elsif ( $exp == 5 )  { return 1024; }
    elsif ( $exp == 6 )  { return 4096; }
    elsif ( $exp == 7 )  { return 16384; }
    elsif ( $exp == 8 )  { return 65536; }
    elsif ( $exp == 9 )  { return 262144; }
    elsif ( $exp == 10 ) { return 1048576; }
  }
  my $res = $num;
  for ( my $i = 2 ; $i <= $exp ; $i++ ) { $res *= $num; }
  return $res;
}

=head3 log_k

Function for calculate the logarithm of any base.

Call: log_k( $base, $num ) NUMBER, NUMBER

Return: $res NUMBER

=cut

sub log_k
{
  my $base = shift @_;
  my $num  = shift @_;
  my $res  = 0;
  if ( $num > 0 )
  {
    $res = log($num) / log($base);
  } else
  {
    die "Cannot calculate log_k($base, $num)\n";
  }
  return $res;
}

=head3 countWords

Function for count words in a sequence.

Call: countWords( $seq, $word ) STRING, NUMBER

Return: %results HASH (KEYS are the elements)

=cut

sub countWords
{
  my $seq   = shift @_;
  my $word  = shift @_;
  my $len   = length $seq;
  my %count = ();

  # Init state when word == 1
  foreach my $b ( @alphabet )
  {
    $count{$b} = 0;
  }

  for ( my $i = 0 ; $i <= ( $len - $word ) ; $i++ )
  {
    my $elem = substr( $seq, $i, $word );
    next if ( $elem =~ /[^ACGT]/ );
    $count{$elem}++;
  }
  return %count;
}

=head3 createWords

Function to create a list of possible words with a specific length

Call: createWords( $kmer, @alphabet ) NUMBER, ARRAY

Return: @old ARRAY

=cut

sub createWords
{
  my $kmer = shift @_;
  $kmer--;
  my @old = @_;
  my @new = ();
  if ( $kmer < 1 )
  {
    return @old;
  } else
  {
    foreach my $e ( @old )
    {
      foreach my $n ( @alphabet )
      {
        push @new, "$e$n";    # add new element
      }
    }
    createWords( $kmer, @new );    # recursion call
  }
}

=head3 randSeq

Function to create a random DNA sequence

Call: randSeq( $size ) NUMBER

Return: $seq STRING

=cut

sub randSeq
{
  my $size = shift @_;
  $size ||= 1000;
  my $seq = '';
  for ( my $i = 0 ; $i <= $size ; $i++ )
  {
    $seq .= $alphabet[ int( rand @alphabet ) ];
  }
  return $seq;
}
1;

=head1 SEE ALSO


=head1 AUTHOR

Juan Caballero, E<lt>jcaballero@systemsbiology.orgE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2009 by Juan Caballero

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.8 or,
at your option, any later version of Perl 5 you may have available.


=cut

__END__
