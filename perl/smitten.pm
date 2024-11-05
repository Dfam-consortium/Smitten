#!/usr/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) smitten.pm
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      This is a reference implementation of the Smitten sequence 
##      identifier format.
##
#******************************************************************************
#*  This software is provided ``AS IS'' and any express or implied            *
#*  warranties, including, but not limited to, the implied warranties of      *
#*  merchantability and fitness for a particular purpose, are disclaimed.     *
#*  In no event shall the authors or the Institute for Systems Biology        *
#*  liable for any direct, indirect, incidental, special, exemplary, or       *
#*  consequential damages (including, but not limited to, procurement of      *
#*  substitute goods or services; loss of use, data, or profits; or           *
#*  business interruption) however caused and on any theory of liability,     *
#*  whether in contract, strict liability, or tort (including negligence      *
#*  or otherwise) arising in any way out of the use of this software, even    *
#*  if advised of the possibility of such damage.                             *
#*                                                                            *
#******************************************************************************
#
=head1 NAME

smitten - Reference implementation of the Smitten sequence identifier format

=head1 SYNOPSIS

  use smitten;

  my ( $assemblyID, $sequenceID, $rangesRef ) = 
             smitten::parseID( "hg38:chr1:100-200_+", 0 );

  my $nID = smitten::normalizeID( "hg38:chr1:100-200_-:10-20_+" );

=head1 DESCRIPTION

Library of functions to manipulate sequence identifiers in the Smitten 
format.  This includes legacy versions of the formats used by Arian
Smit in his various tools.
   
For example to normalize a sequence identifier in V0, V1, or
V2 Smitten format:

Smitten V0:  
      chr1_100_200          positive strand range
      hg38:chr1_100_200     positive strand range with assembly ID
      hg38:chr1_100_200_R   negative strand range

Smitten V1:
      chr1:1-200            forward strand range, seq identifier
      chr1:200-1            reverse strand range, seq identifier

Smitten V2:
      chr1:100-200_+        positive strand range
      hg38:chr1:100-200_+   positive strand range with assembly ID
      hg38:chr1:100-200_-   negative strand range 

In addition, ranges may be chained. For example:

      hg38:chr1:100-200_+:10-50_-:1-5_+   

is a V2 format that is equivalent to:
    
      hg38:chr1:145-149_-

in its normalized form. To parse these various formats reliably this library
includes a parser that can recognize all three.  For example:

     my ( $assemblyID, $sequenceID, $rangesRef ) = 
             smitten::parseID( "hg38:chr1:100-200_+", 0 );

will return:
     $assemblyID = "hg38";
     $sequenceID = "chr1";
     $rangesRef  = A reference to a list of start, end, orientation tuples:
                     [ 
                       [ 10, 200, '+'] 
                     ]

The chained identifier:
     my ( $assemblyID, $sequenceID, $rangesRef ) = 
             smitten::parseID( "chr1:100-200_+:10-30_-" );

would generate the following ranges structure:
                     [ 
                       [ 10, 30, '-'] 
                       [ 100, 200, '+' ],
                     ]

Finally, chained identifiers may also be normalized with the normalizeID() function:

   my $nID = smitten::normalizeID( "hg38:chr1:100-200_-:10-20_+" );

which represents the positions 10 to 20 in the reversed sequence from
100 to 200.  This routine would return "hg38:chr1:180-190_-".

   
The options are:


=over 4

=item -version

Displays the version of the program

=back

=head1 SEE ALSO

=head1 COPYRIGHT

Copyright 2022-2024 Robert Hubley, Institute for Systems Biology

=head1 LICENSE

This code may be used in accordance with the Open Source License v. 3.0
http://opensource.org/licenses/OSL-3.0

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>

=cut

#
# Module Dependence
#
use strict;
use Getopt::Long;
use Data::Dumper;
use Carp;
use FindBin;

my $Version = "0.1";
my $DEBUG = 0;

##-------------------------------------------------------------------------##

=head2 convertID()

  Use:  my $v2SeqID = 
             smitten::convertID( id => 'chr1:1-10_+',
                                [zbho => 0/1],
                               );

      id : Sequence identifier in V0, V1, or V2 Smitten format.

     zbho: For conversion of legacy identifiers, this treats
               the coordinate ranges as Zero-Based, Half-Open
               coordinates.
  
  Ranges are optional and may be in V0, V1 or V2 formats. Beware
  that these formats are ambiguous in some cases.  For instance,
  a sequence ID of "chr1_1_1" is allowed as a non-ranged sequence ID in
  V2 format but is indistinguishable from a sequence ID of "chr1" from 
  position 1 to 1 in V0 format. This routine infers the format by 
  identifying sequence range suffixes from the identifier string in either
  the V0, V1, or V2 syntax.  The matching range is removed from the identifier
  and the remaining identifier is recursively searched until either no range
  pattern is found, or version mismatch is identified.  In the above example,
  the sequence ID "chr1_1_1" would be identified as a V0 format sequence ID,
  and the V2 equivalent would be returned as "chr1:1-1_+".

  Returns:

       seqID   : The equivalent Smitten V2 format sequence 
                 identifier.
      idVersion: The inferred version of the sequence identifier:
                   0: V0 format
                   1: V1 format
                   2: V2 format

=cut

##-------------------------------------------------------------------------##·
sub convertID {
  my (%args) = ( 'zbho' => 0,
                 @_    
               );

  croak "convertID: missing 'id' argument!\n"
      unless ( exists $args{'id'}  );

  my $idStr = $args{'id'};
  my $isZeroBasedHalfOpen = $args{'zbho'};
  my $parseFmt = $args{'format'};
  my $origIDStr = $idStr;
  my $idPrefix = $idStr;
  my @ranges = ();

  croak "convertID: Identifier \'$origIDStr\' contains a space or a line termination character!\n"
    if ( $origIDStr =~ /[\s\n\r]/ );

  my $inferredFmt = "";
  # Process ranges from right to left 
  #   E.g. chr1:1-100:5-10_+
  #            loop 1: ":5-10_+"
  #            loop 2: ":1-100"
  while ( $idStr =~ s/(.*)(([:_])(\d+)([-_])(\d+)((_)([R\+\-]))?)$/$1/ )
  {
    my $rangeStr = $2;
    my $start = $4;
    my $end = $6;
    my $orient = $9;

    # Infer the format version of the range
    my $rangeFmt;
    if ( $3.$5.$8 eq ":-_" ) {
      $rangeFmt = 2;
    }elsif ( $3.$5 eq "__" ) {
      $rangeFmt = 0;
    }elsif ( $3.$5 eq ":-" ) { 
      $rangeFmt = 1;
    }else {
      $idStr .= $rangeStr;
      last;
    }

    # Catch mismatch in range formats or when strict
    # conversion is used.
    #    E.g.  chr1_1_100:5-10_+
    #          This can be either:
    #             Format     SeqID               Ranges
    #               V2      chr1_1_100           [5,10,'+']
    #               V1      chr1_1_100:5-10_+
    #               V0      chr1_1_100:5-10_+
    #          V2 is chosen based on the suffix identification
    #          and the remaining incompatable string remains
    #          as the inferred sequence identifier.
    if ( $inferredFmt ne "" && $rangeFmt != $inferredFmt ) {
       # Place mismatched range back on the sequence identifier
       $idStr .= $rangeStr;
       last;
    } 
    $inferredFmt = $rangeFmt;
    
    if ( $start > $end ) {
      if ( $inferredFmt != 1 ) {
        croak "convertID: V$inferredFmt identifier \'$origIDStr\' must have increasing range order!\n";
      }
      $orient = '-';
      my $tmp = $end;
      $end = $start;
      $start = $tmp;
    }elsif ( $orient eq "" ) {
      $orient = '+';
    }
    $orient = '-' if ( $orient eq 'R' );

    # Convert back into one-based, fully closed
    $start++ if ( $isZeroBasedHalfOpen );

    if ( $start < 1 || $end < 1 ) {
      croak "convertID: The converted range ($start-$end) contains coordinate that is less " . 
            "than one. ( id=\'$origIDStr\', isZeroBasedHalfOpen=$isZeroBasedHalfOpen ).\n";
    }

    unshift @ranges, [$start, $end, $orient];
  }
  # validate subranges
  for ( my $i = 1; $i <= $#ranges; $i++ ) {
    my $prevLen = $ranges[$i-1]->[1] - $ranges[$i-1]->[0] + 1;
    croak "convertID: Sequence sub-range outside parent range ( \'$origIDStr\' )."
        if  ($ranges[$i]->[0] > $prevLen || $ranges[$i]->[1] > $prevLen );
  }

  # validate the assembly/sequence identifier(s)
  # ":" returns empty @ids
  my @ids = split(/:/, $idStr);
  my $colonCount = $idStr =~ tr/:/:/;

  my $assemblyID = "";
  my $sequenceID = "";
  if ( $colonCount == 1 && @ids == 2 && $ids[0] ne "" && $ids[1] ne "") {
    $assemblyID = $ids[0];
    $sequenceID = $ids[1];
  }elsif ( $colonCount == 0 && $idStr ne "" ) { 
    $sequenceID = $idStr;
  }else {
    # More than one ':' in the identifier -- cannot create a valid V2 out of this
    croak "convertID: Identifier \'$origIDStr\' contains more than one ':' in it's assembly+sequence identifiers.\n";
  }
  
  if ( $sequenceID eq "" ) {
    croak "convertID: Identifier \'$origIDStr\' does not have a sequence identifier!\n";
  }

  my $v2ID;
  $v2ID .= $assemblyID.":" if ( $assemblyID ne "" );
  $v2ID .= $sequenceID;
  foreach my $range ( @ranges ) {
    $v2ID .= ":" . $range->[0] . "-" . $range->[1] . "_" . $range->[2];
  }
  
  return ($v2ID, $inferredFmt);
}


##-------------------------------------------------------------------------##

=head2 parseID()

  Use:  my ( $assemblyID, $sequenceID, $rangesRef ) = 
                                     smitten::parseID( $idStr );

      $idStr : Sequence identifier V2 Smitten format

  Ranges are optional but must be in V2 format if present.

  Returns the components of a Smitten identifier as:

      $assemblyID : The assembly component as a string or an empty string
                    if not defined.
      $sequenceID : The sequence identifier as a string.
      $rangesRef  : A reference to a list of start, end, orientation tuples.
                    The first entry in the list represents the deepest range.
                    E.g: chr1:100-200_+:10-30_- would generate the following
                    ranges structure:

                           [ 
                             [ 10, 30, '-'] 
                             [ 100, 200, '+' ],
                           ]

=cut

##-------------------------------------------------------------------------##·
sub parseID {
  my $idStr = shift;

  my $origIDStr = $idStr;
  my @ranges = ();

  croak "parseID: identifier \'$origIDStr\' contains a space or a line termination character!\n"
    if ( $idStr =~ /\s\n\r/ );

  while ( $idStr =~ s/(.*)(([:])(\d+)([-])(\d+)((_)([\+\-]))?)$/$1/ )
  {
    my $rangeStr = $2;
    my $start = $4;
    my $end = $6;
    my $orient = $9;
    
    if ( $start > $end ) {
      croak "parseID: V2 identifiers \'$origIDStr\' must have increasing range order!\n";
    }elsif ( $orient eq "" ) {
      croak "parseID: V2 identifiers \'$origIDStr\' must have orientation specified for every subrange!\n";
    }
    unshift @ranges, [$start, $end, $orient];
  }
  # validate subranges
  for ( my $i = 1; $i <= $#ranges; $i++ ) {
    my $prevLen = $ranges[$i-1]->[1] - $ranges[$i-1]->[0] + 1;
    croak "parseID: Sequence sub-range outside parent range ( $origIDStr )."
        if  ($ranges[$i]->[0] > $prevLen || $ranges[$i]->[1] > $prevLen );
  }

  my @ids = split(/:/, $idStr);

  my $assemblyID = "";
  my $sequenceID = "";
  if ( @ids == 2 ) {
    $assemblyID = $ids[0];
    $sequenceID = $ids[1];
  }elsif ( @ids == 1 ) {
    $sequenceID = $ids[0];
  }
  
  if ( $sequenceID eq "" ) {
    croak "parseID: Identifier \'$origIDStr\' does not have a sequence identifier!\n";
  }
  
  return( $assemblyID, $sequenceID, \@ranges );
}

##-------------------------------------------------------------------------##

=head2 normalizeID()

  Use:  my $normalizedIDStr = smitten::normalizeID( $idStr );

      $idStr : Sequence identifier in V2 Smitten format.

  Returns the normalized smitten identifier consisting of a single
  sequence range given an identifier that has more than one level of 
  sequence ranges.  For example, given the input identifier:

          hg38:chr1:100-200_-:10-20_+

  representing the positions 10 to 20 in the reversed sequence from
  100 to 200.  This routine would return:

          hg38:chr1:180-190_-

=cut

##-------------------------------------------------------------------------##·
sub normalizeID {
  my $idStr = shift;

  my ( $assemblyID, $sequenceID, $rangesRef ) = parseID( $idStr );
  
  # Nothing to normalize
  return $idStr if ( @{$rangesRef} == 0 );
 
  my $startIdx = $rangesRef->[$#{$rangesRef}]->[0];
  my $endIdx = $rangesRef->[$#{$rangesRef}]->[1];
  my $currOrient = $rangesRef->[$#{$rangesRef}]->[2];
  for ( my $i = $#{$rangesRef}-1; $i >= 0; $i-- ) {
    if ( $rangesRef->[$i]->[2] eq "-" ) {
      $startIdx = $rangesRef->[$i]->[1]-$startIdx+1;
      $endIdx = $rangesRef->[$i]->[1]-$endIdx+1;
    }else {
      $startIdx = $rangesRef->[$i]->[0]+$startIdx-1;
      $endIdx = $rangesRef->[$i]->[0]+$endIdx-1;
    }
    if ( $rangesRef->[$i]->[2] eq '-' && $currOrient eq '-' ) {
      $currOrient = '+';
    }elsif ( $rangesRef->[$i]->[2] eq '-' || $currOrient eq '+' ) {
      $currOrient = '-';
    }
  }

  my $retStr = "";
  $retStr .= $assemblyID . ":" if ( $assemblyID ne "" );
  $retStr .= $sequenceID . ":";
  if ( $startIdx < $endIdx ) {
    $retStr .= $startIdx . "-" . $endIdx;
  }else {
    $retStr .= $endIdx . "-" . $startIdx;
  }
  $retStr .= "_$currOrient";

  return $retStr;

}
 

##-------------------------------------------------------------------------##
## Use: my _privateMethod( $parameter => value );
##
##      $parameter       : A parameter to the method
##
##  Returns
##      Methods with the prefix "_" are conventionally considered
##      private.  This bit-o-documentation is not formatted to
##      print out when perldoc is run on this file.
##
##-------------------------------------------------------------------------##
sub _privateMethod {
  my %parameters = @_;

  my $subroutine = ( caller( 0 ) )[ 0 ] . "::" . ( caller( 0 ) )[ 3 ];
  print "$subroutine( " .  @{[%parameters]} . "): Called\n" if ( $DEBUG );

}


##-------------------------------------------------------------------------##

=head2 publicMethod()

  Use:  my $retVal = publicMethod( $parameter1 => value, 
                                   $parameter2 => value );

    $parameter1:   A generic scalar parameter
    $parameter2:   A generic scalar parameter

  $retVal contains the scalar result of this subroutine.  This
  is a public function and this documentation will print out
  when perldoc is run on this file.

=cut

##-------------------------------------------------------------------------##·
sub publicMethod {
  my %parameters = @_;

  my $subroutine = ( caller( 0 ) )[ 0 ] . "::" . ( caller( 0 ) )[ 3 ];
  print "$subroutine( " .  @{[%parameters]} . "): Called\n" if ( $DEBUG );
}

1;
