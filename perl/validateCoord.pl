#!/usr/local/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) validateCoord.pl
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      This is a template for generic perl scripts.  It
##      includes an area for POD documentation (ie. perldoc this-file )
##      and a generic way to store startup parameters in a file
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

validateCoord.pl - validate Smitten coordinate sequences against a reference

=head1 SYNOPSIS

  validateCoord.pl [-version] [-repair|r <repaired file>] 
                   -twobit <twoBit reference file> 
                   <stk or fasta file>

=head1 DESCRIPTION

Validate the Smitten coordinates stored in a FASTA, or Stockholm file.
Optionally repair the coordinates if they are fixable.

TODO: Expand to CSV/TSV with field designator.  Perhaps also use the 
sequence validator to validate other identifier formats as well -- i'm
looking at you BED/GFF3.

The options are:

=over 4

=item -version

Displays the version of the program

=back

=head1 SEE ALSO

=head1 COPYRIGHT

Copyright 2022 Robert Hubley, Institute for Systems Biology

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
use FindBin;
use Carp;
use File::Temp qw/ tempfile tempdir /;
use lib $FindBin::RealBin;
use smitten;

#
# If libraries are to be bundled use FindBin to
# locate and load the libraries.
#
# e.g. for a module in the same directory the following should resolve it:
#   use lib "$FindBin::RealBin";
#   use MyLibInSameDir;
#
# or if it's in a directory above:
#
#   use lib "$FindBin::RealBin/..";
#   use MyLibInDirAbove
#

my $Version = "0.1";

# Magic numbers/constants here
my $DEBUG = 0;
#
# The path to the dependencies 'twoBitInfo' and
# 'twoBitToFa' that are part of the UCSC Tools package.
#
# TODO: Formalize setting this location for UCSC Tools
my $pathToUCSCTools = "/usr/local/ucscTools";
#
my $pathToCopmem2 = "/usr/local/src/copmem2-2.0";

#
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number paramters
#
my @getopt_args = (
    '-version', # print out the version and exit
    '-twobit|t=s',
    '-fasta|f=s',
    '-copmem2|c',
    '-repair|r=s',
);

my %options = ();
Getopt::Long::config("noignorecase", "bundling_override");
unless (GetOptions(\%options, @getopt_args)) {
    usage();
}

sub usage {
  print "$0 - $Version\n\n";
  exec "pod2text $0";
  exit;
}

if ($options{'version'}) {
  print "$Version\n";
  exit;
}

#
# ARGV Processing
#
if ( !@ARGV  ) {
  usage();
}

if ( !defined $options{'twobit'} && !defined $options{'fasta'}) {
  die "You must specify a reference 2bit file with the -twobit option!\n";
}
my $twoBit = $options{'twobit'};

my $seqFile = $ARGV[0];

open IN,"<$seqFile" or die "Could not open $seqFile";
my $maxLines = 100;
my $fType = "Unknown";
while (<IN>){
  last if ( $maxLines-- == 0 );
  if ( /^#\s*STOCKHOLM\s+\d/ ) {
    $fType = "stk";
    last;
  }
  if ( /^>\S+/ ) {
    $fType = "fasta";
    last;
  }
}
close IN;
if ( $fType eq "Unknown" ) {
  die "Could not identify input file type.  Should be either fasta or stockholm.\n";
}

my @seqRecs = ();
my $IN;
open $IN, "<$seqFile" or die "Could not open $seqFile for reading!\n";
if ( $fType eq "stk" ) {
  my $familyID;
  while (<$IN>){
    # TODO: Optionally add familyID/Accession to sequence records, for improved error reporting
    $familyID = $1 if ( /^#=GF\s+ID\s+(\S+)/ );
    next if ( /^#/ || /^\s*$/ || /^\/\// );
    if ( /(\S+)\s+([ACGTRYWSKMNXBDHV\.\-]+)\s*$/ ) {
      my $id = $1;
      my $seq = $2;
      $seq =~ s/[\.\-]//g;
      push @seqRecs,[$id, uc($seq)];
    }
  }
  close $IN;
}elsif ( $fType eq "fasta" ) {
  my $seq;
  my $id;
  while (<$IN>){
    if ( />(\S+)/ ) {
      my $tmpID = $1;
      if ( $seq ) {
          push @seqRecs,[$id, uc($seq)];
      }
      $seq = "";
      $id = $tmpID;
      next;
    }
    s/[\n\r\s]+//g;
    $seq .= $_;
  }
  if ( $seq ) {
    push @seqRecs,[$id, uc($seq)];
  }
}
close $IN;

my $validationErrors;
if ( $options{'fasta'} ) {
  $validationErrors = &_validateSequencesCopmem2( $options{'fasta'}, \@seqRecs, $pathToCopmem2 );
}else {
  $validationErrors = &_validateSequences( $twoBit, \@seqRecs, $pathToUCSCTools );
}

print "$seqFile: Found " . scalar(@{$validationErrors}) . " errors in sequence set!\n";
if ( @{$validationErrors} ) {
  &_interpretDiffs($validationErrors);

  # TODO Should we generate a repaired file if there are no errors...currently no
  if ( $options{'repair'} ) {
    my $repairFile = $options{'repair'};
    my $OUT;
    my $IN;
    open $OUT,">$repairFile" or die "Could not open $repairFile for writing!\n";
    open $IN,"<$seqFile" or die "Could not open $seqFile";
    my %errorIndices = ();
    my $errorIdx = 0;
    foreach my $error ( @{$validationErrors} ) {
      $errorIndices{$error->[0]} = $errorIdx; 
      $errorIdx++;
    }
    if ( $fType eq "stk" ) {
      my $seqIdx = 0;
      while ( <$IN> )  {
        if ( /^#/ || /^\s*$/ || /^\/\// ) {
          print $OUT $_;
          next;
        }
        if ( /(\S+)\s+([ACGTRYWSKMNXBDHV\.\-]+)\s*$/ ) {
          my $id = $1;
          my $seq = $2;
          if ( exists $errorIndices{$seqIdx} ) {
            my $fixedID = &fixID( $id, $validationErrors->[$errorIndices{$seqIdx}] );
            print $OUT "$fixedID $seq\n";
          }else {
            my $fixedID = &fixID( $id );
            print $OUT "$fixedID $seq\n";
          }
          $seqIdx++;
        }
      }
    }else {
      die "NOT IMPLEMENTED YET!\n";
    }
    close $IN;
    close $OUT;
  }
}

sub fixID {
  my $origID = shift;
  my $error = shift;

  # Assuming 1-based, fully closed ( provide an option to override that )
  my ($id, $version);
  eval {( $id, $version ) = convertID( id => $origID );};
  if ( $@ ) {
    eval {( $id, $version ) = convertID( id => $origID, zbho=>1 );};
    if ( $@ ) {
      print "$File::Find::name : $origID: produced an error $@\n";
    }
  }
  my $normID = normalizeID($id);
  my ($assembly, $seqID, $ranges) = parseID($normID,0);
  my $start = $ranges->[0]->[0];
  my $end = $ranges->[0]->[1];
  my $orient = $ranges->[0]->[2];


  if ( $error ) {
    $start += $error->[2];
    $end += $error->[3];
    if ( $error->[4] ne "" ) {
      $orient = $error->[4];
    }
  }

  my $newID;
  $newID = $assembly . ":" if ( $assembly );
  $newID = $seqID . ":" . $start . "-" . $end . "_" . $orient;
 
  return $newID;
}


######################## S U B R O U T I N E S ############################


##-------------------------------------------------------------------------##
## 
## _interpretDiffs()
##
##-------------------------------------------------------------------------##
sub _interpretDiffs {
  my $validationErrors = shift;

  foreach my $error ( @{$validationErrors} ) {
     my $seqIdx = $error->[0];
     my $origID = $error->[1];
     my $startOffset = $error->[2];
     my $endOffset = $error->[3];
     my $orientationChange = $error->[4];
     my $altMatches = $error->[5];
     print "[$seqIdx] $origID => ";

     # Assuming 1-based, fully closed ( provide an option to override that )
     my ($convID, $version);
     eval {( $convID, $version ) = convertID( id => $origID );};
     if ( $@ ) {
       eval {( $convID, $version ) = convertID( id => $origID, zbho=>1 );};
       if ( $@ ) {
         print "$File::Find::name : $origID: produced an error $@\n";
       }
     }
     my $normID = normalizeID($convID);
     my ($assembly, $seqID, $ranges) = parseID($normID,0);
     my $oldStart = $ranges->[0]->[0];
     my $oldEnd = $ranges->[0]->[1];
     my $orient = $ranges->[0]->[2];

     if ( $orientationChange ne "" ) {
       $orientationChange = ", and orientation changed to \'$orientationChange\'";
       $orient = $orientationChange;
     }
     if ( $startOffset == 1 && $endOffset == -1 ) {
       print "" . ( $oldStart + $startOffset ) . "-" . ( $oldEnd + $endOffset ) ."_". $orient . " ";
       print " Invalid coordinate system - zeroBasedHalfOpen$orientationChange\n";
     }elsif ( $startOffset == 1 && $endOffset == 0 ) {
       print "" . ( $oldStart + $startOffset ) . "-" . ( $oldEnd + $endOffset ) ."_". $orient . " ";
       print " Invalid coordinate system - zeroBasedFullyClosed$orientationChange\n";
     }elsif ($startOffset == 0 && $endOffset == -1 ) {
       print "" . ( $oldStart + $startOffset ) . "-" . ( $oldEnd + $endOffset ) ."_". $orient . " ";
       print " Invalid coordinate system - oneBasedHalfOpen$orientationChange\n";
     }elsif ($startOffset == 0 && $endOffset == 0 && $orientationChange ne "" ) {
       print "" . ( $oldStart + $startOffset ) . "-" . ( $oldEnd + $endOffset ) ."_". $orient . " ";
       print " Invalid orientation!\n";
     }elsif ( $startOffset ne "" ) {
       print "" . ( $oldStart + $startOffset ) . "-" . ( $oldEnd + $endOffset ) ."_". $orient . " ";
       print " Invalid coordinate system - shifted by start:$startOffset, end:$endOffset" 
             . $orientationChange . "\n";
     }elsif ( @{$altMatches} ) {
       print "" . ( $oldStart + $startOffset ) . "-" . ( $oldEnd + $endOffset ) ."_". $orient . " ";
       print " Alt matches! [offset,orientation]: ";
       my @strs = ();
       foreach my $match ( @{$altMatches} ) {
          push @strs,"[" . join(",", @{$match} ) . "]";
       }
        print "" . join(",",@strs) . "\n";
     }else {
       print " Could not find a match!\n";
     }
  }
}


##-------------------------------------------------------------------------##
##
## _validateSequences()
##
##  Use:  my $validationErrors = validateSequences( $twoBitRef
##                            $seqRecs, $pathToUCSCTools, [$maxPadding] );
##
##    $twoBitRef:       Path to twobit sequence file containing the
##                        reference sequences for the validation.
##    $seqRecs:         A reference to a list of list references containing
##                        the sequences and idenfiers that are to be validated.
##                        See below for structure.
##    $pathToUCSCTools: A path to the directoy containing the UCSC tools 
##                        'twoBitInfo' and 'twoBitToFa'
##    $maxPadding:      Optional amount (bp) of flanking sequence to use
##                        while validating the coordinates.  (default:5)
##
##  This routine validates coordinate-sequence pairs where the coordinates
##  have been encoded in the smitten format.  The comparison is made to a
##  reference sequence file in the UCSC 2bit format.  Sequence identifier
##  pairs are passed in as a list of lists reference in the form:
##
##     [ 
##        [ id1, seq1 ],
##        [ id2, seq2 ],
##        [ id3, seq3 ],
##        ...
##     ]
##
##  And errors obtained from the validation are returned as a list of 
##  records in the form:
##
##    [
##       [ seqIdx, origID, startOffset, endOffset, 
##         orientChange, [ [startOffset1, endOffset1, oChange1], ...]],
##       [ seqIdx, origID, startOffset, endOffset, 
##         orientChange, [ [startOffset1, endOffset1, oChange1], ...]],
##       ...
##    ]
##  
##  Where:
##
##      seqIdx:       The zero-based index into seqRecs of the sequence 
##                      failing validation.
##      origID:       The original smitten sequence identifier failing
##                      validation.
##      startOffset:  The offset from the original start position for
##                      a sequence with a single match.
##      endOffset:    The offset from the original end position for a
##                      sequence with a single match.
##      orientChange: The new orientation for the sequence, if it changed
##                      ('-' or '+') or an empty string if it didn't change.
##      altMatches:   A reference to a list of alternative matching locations
##                      if a single match wasn't found.  The alternative matching
##                      locations are stored as lists of startOffset, endOffset
##                      and orientChange fields as previously described.
##
##-------------------------------------------------------------------------##
sub _validateSequences {
  my $twoBitRef = shift;
  my $seqRecs = shift;
  my $pathToUCSCTools = shift;
  my $maxPadding = shift;

  $maxPadding = 5 if ( ! defined $maxPadding || $maxPadding eq "" || $maxPadding < 0 );

  # Open up the two bit file and get seq sizes
  open IN,"$pathToUCSCTools/twoBitInfo $twoBitRef stdout|" or 
     croak "validateSequences: Error could not open $twoBitRef with $pathToUCSCTools/twoBitInfo!\n";
  my %seqSize = ();
  while ( <IN> ) {
    if ( /^(\S+)\t(\d+)\s*$/ ) {
      $seqSize{$1} = $2;
    }
  }
  close IN;

  # Generate a BED file containing the sequence ranges we wish
  # to lookup each padded so that we can tolerate errors in the
  # sequence coordinates.
  my ( $tmpFH, $tmpFilename ) =
    tempfile( UNLINK => 0, SUFFIX => ".bed", DIR => "." );
  my @leftPads = ();
  my @rightPads = ();
  my @seqRanges = ();
  my $idx = 0;
  foreach my $seqRec ( @{$seqRecs} ){
    my $origID = $seqRec->[0];
    # Assuming 1-based, fully closed ( provide an option to override that )
    my ($convID, $version);
    eval {( $convID, $version ) = convertID( id => $origID );};
    if ( $@ ) {
      eval {( $convID, $version ) = convertID( id => $origID, zbho=>1 );};
      if ( $@ ) {
        print "$File::Find::name : $origID: produced an error $@\n";
      }
    }
    my $normID = normalizeID($convID);
    my ($assembly, $seqID, $ranges) = parseID($normID,0);
    push @seqRanges, [$seqID, $ranges->[0]->[0], $ranges->[0]->[1], $ranges->[0]->[2]];
    # Convert to 0-based, half open and allow for variable pads
    my $adjustedStart = $ranges->[0]->[0]-1; 
    if ( $adjustedStart >= $maxPadding ) {
      $adjustedStart -= $maxPadding;
      push @leftPads,$maxPadding;
    }else {
      push @leftPads,$adjustedStart;
      $adjustedStart = 0;
    }
    # No drecrement as this is an open coordinate
    my $adjustedEnd = $ranges->[0]->[1];
    if ( $adjustedEnd <= $seqSize{$seqID} - $maxPadding ) {
      $adjustedEnd += $maxPadding;
      push @rightPads, $maxPadding;
    }else {
      push @rightPads,($seqSize{$seqID} - $adjustedEnd);
      $adjustedEnd = $seqSize{$seqID};
    }
    print $tmpFH "$seqID\t$adjustedStart\t$adjustedEnd\tseq-$idx\n";
    $idx++;
  }
  close $tmpFH;

  # Extract sequences from the reference 2bit file
  open IN,"$pathToUCSCTools/twoBitToFa -bed=$tmpFilename $twoBitRef stdout|" or 
     croak "validateSequences: Error could not open $twoBitRef with $pathToUCSCTools/twoBitToFa!\n";
  my $tmp_seq;
  my $seq = "";
  my @refSeqs = ();
  $idx = "";
  while (<IN>){
    if ( />seq-(\d+)/ ) {
      my $tIdx = $1;
      $refSeqs[$idx] = uc($seq) if ( $seq ne "" );
      $seq = "";
      $idx = $tIdx;
      next;
    }
    s/[\n\r\s]+//g;
    $seq .= $_;
  }
  $refSeqs[$idx] = uc($seq) if ( $seq ne "" );
  close IN;
  unlink($tmpFilename);

  # Validate sequences/ranges against the reference set
  my @validationErrors = ();
  for ( my $i = 0; $i <= $#seqRanges; $i++ ){
    # Unpack datastructures
    my $idRec = $seqRanges[$i];  
    my $start = $idRec->[1];
    my $end = $idRec->[2];
    my $orient = $idRec->[3];
    # The length of the pads acheived on this particular 
    # sequence. This may be less than the padding specified
    # if the at the begining or end of a contig.
    my $leftPad = $leftPads[$i];
    my $rightPad = $rightPads[$i];
    # The original sequence and un-parsed identifier
    my $origID = $seqRecs->[$i]->[0];
    my $origSeq = $seqRecs->[$i]->[1];
    # and its reverse complement
    my $revOrigSeq = reverse($origSeq);
    $revOrigSeq =~ tr/ACGTRYWSKMNXBDHV/TGCAYRSWMKNXVHDB/;

    # The reference sequence pulled from the twobit file
    my $refSeq = $refSeqs[$i];

    # Identify matches by shifting over all possible
    # origSeq-length substrings of the padded reference
    # sequence.
    my $lenSeq = length($origSeq);
    my $lenDiff = length($refSeq) - $lenSeq;
    my @possibleMatches = ();
    my $startOffset = '';
    my $endOffset = '';
    my $orientStr = '';
    my $validMatch = 0;
    for ( my $j = 0; $j <= $lenDiff; $j++ ) {
      # Convert this $j index into offsets from the
      # original start and end sequence indexes.
      my $sOffset = $j - $leftPad;
      my $eOffset = $sOffset + ($lenSeq-($end-$start+1));

      # Compare subsequences in both orientations
      my $matchDir = '';
      if ( $orient eq "+" ) {
        # Assume '+' first
        if ( $origSeq eq substr($refSeq,$j,$lenSeq) ) {
          $matchDir = '+'; # Match to positive strand
        }elsif ( $revOrigSeq eq substr($refSeq,$j,$lenSeq) ) {
          $matchDir = '-'; # Match to reverse strand
        }
      }else {
        # Assume '-' first
        if ( $revOrigSeq eq substr($refSeq,$j,$lenSeq) ) {
          $matchDir = '-'; # Match to reverse strand
        }elsif ( $origSeq eq substr($refSeq,$j,$lenSeq) ) {
          $matchDir = '-'; # Match to positive strand
        }
      }

      # If we have a match
      if ( $matchDir ne '' ) {
        # If the match is where we were expecting to find it
        # (start offset == 0), process this as the preferred 
        # hit.
        if ( $sOffset == 0 ) {
          # Correct-ish.  It could be in the reverse, 
          # orientation, or at a different end offset.
          $startOffset = $sOffset;
          $endOffset = $eOffset;
          if ( $matchDir ne $orient ) {
            if ( $orient eq "+" ) {
              $orientStr = "-";
            }else {
              $orientStr = "+";
            }
          }else {
            # This flags that this was a perfect match
            $validMatch = 1 if ( $eOffset == 0 );
          }
        }else {
          # If there was a match and it's not where we expected
          # (start offset == 0), put this on a list of alternative
          # matches.
          push @possibleMatches, [$sOffset,$eOffset,$matchDir];
        }
      }
    }
    next if ( $validMatch );

    if ( $startOffset ne '' ) {
      push @validationErrors, [$i, $origID, $startOffset, $endOffset, $orientStr, []];
    }elsif ( @possibleMatches ) {
      if ( scalar(@possibleMatches) == 1 ) {
        my $orientStr = '';
        if ( $possibleMatches[0]->[2] ne $orient ) {
          $orientStr = $possibleMatches[0]->[2];
        }
        push @validationErrors, [$i, $origID, $possibleMatches[0]->[0], $possibleMatches[0]->[1], $orientStr, []];
      }else { 
        push @validationErrors, [$i, $origID, "", "", "", \@possibleMatches];
      }
    }else {
      push @validationErrors, [$i, $origID, "", "", "", []];
    }
  }

  return ( \@validationErrors );
}

##-------------------------------------------------------------------------##
##
## _validateSequencesCopmem2()
##
##  Use:  my $validationErrors = validateSequencesCopmem2( $fastaRef,
##                            $seqRecs, $pathToCopmem2 );
##
##    $fastaRef:        Path to FASTA sequence file containing the
##                        reference sequences for the validation.
##    $seqRecs:         A reference to a list of list references containing
##                        the sequences and idenfiers that are to be validated.
##                        See below for structure.
##    $pathToCopmem2:   A path to the directoy containing the copmem2 tool.
##
##  This routine validates coordinate-sequence pairs where the coordinates
##  have been encoded in the smitten format.  The comparison is made to a
##  reference sequence file in FASTA format using copmem2 MEM finder.  
##  Sequence identifier pairs are passed in as a list of lists reference 
##  in the form:
##
##     [ 
##        [ id1, seq1 ],
##        [ id2, seq2 ],
##        [ id3, seq3 ],
##        ...
##     ]
##
##  And errors obtained from the validation are returned as a list of 
##  records in the form:
##
##    [
##       [ seqIdx, origID, startOffset, endOffset, 
##         orientChange, [ [startOffset1, endOffset1, oChange1], ...]],
##       [ seqIdx, origID, startOffset, endOffset, 
##         orientChange, [ [startOffset1, endOffset1, oChange1], ...]],
##       ...
##    ]
##  
##  Where:
##
##      seqIdx:       The zero-based index into seqRecs of the sequence 
##                      failing validation.
##      origID:       The original smitten sequence identifier failing
##                      validation.
##      startOffset:  The offset from the original start position for
##                      a sequence with a single match.
##      endOffset:    The offset from the original end position for a
##                      sequence with a single match.
##      orientChange: The new orientation for the sequence, if it changed
##                      ('-' or '+') or an empty string if it didn't change.
##      altMatches:   A reference to a list of alternative matching locations
##                      if a single match wasn't found.  The alternative matching
##                      locations are stored as lists of startOffset, endOffset
##                      and orientChange fields as previously described.
##
##-------------------------------------------------------------------------##
sub _validateSequencesCopmem2 {
  my $fastaRef = shift;
  my $seqRecs = shift;
  my $pathToCopmem2 = shift;
  my $threads = shift;

  # Generate a FASTA file containing the sequences we wish to map.
  my ( $tmpFH, $tmpFilename ) =
    tempfile( UNLINK => 0, SUFFIX => ".fa", DIR => "." );
  my $recIdx = 0;
  my %tooShort = ();
  foreach my $seqRec ( @{$seqRecs} ){
    print $tmpFH ">seq-$recIdx\n" . $seqRec->[1] . "\n";
    if ( length($seqRec->[1]) < 50 ) {
      $tooShort{"seq-".$recIdx} = 1;
    }
    $recIdx++;
  }
  close $tmpFH;

  # Run copmem2 to find the MEMs
  my ( $tmpMEMFH, $tmpMEMFilename ) =
    tempfile( UNLINK => 0, SUFFIX => ".mems", DIR => "." );
  my $threadParam = "";
  $threadParam = "-t $threads" if ( $threads );
  my $cmd = "$pathToCopmem2/copmem2 -l 50 -b $threadParam -o $tmpMEMFilename $fastaRef $tmpFilename 2>1 > /dev/null";
  system($cmd);

  if ( ! -s $tmpMEMFilename ) {
    die "Error: copmem2 failed to produce any output!\n";
  }

  open IN,"<$tmpMEMFilename" or die "Could not open $tmpMEMFilename for reading!\n";
  my $tmpID;
  my $seqRecsIdx;
  my $orient;
  my @matches = ();
  while (<IN>) {
    #> JAHSPW020000001.1:31258064-31256478
    #> JAHSPW020000001.1:31258064-31256478 Reverse
    if ( /^>\s*(\S+)\s+(Reverse)?/ ) {
      $tmpID = $1;
      if ( $tmpID =~ /seq-(\d+)/ ) {
        $seqRecsIdx = $1;
      }else {
        die "Error: Could not parse sequence index from $tmpID\n";
      }
      if ( $2 ) {
        $orient = "-";
      }else {
        $orient = "+";
      }
      next;
    }
    # ENA|JAHSPW020000001|JAHSPW020000001.1	30375752	26	57
    if ( /^\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s*$/ ) {
      my $genSeqID = $1;
      my $genStart = $2;
      my $queryStart = $3;
      my $matchLen = $4;

      my $seqRec = $seqRecs->[$seqRecsIdx];
      my $seqRecID = $seqRec->[0];
      my $seqRecLen = length($seqRec->[1]);

      if ( $tooShort{$tmpID} ) {
        if ( $orient eq "+" ) {
          # only print once
          print "[$seqRecID] Sequence too short to validate\n";
        }
      }else {
        if ( $matchLen == $seqRecLen ) {
          print "[$seqRecID] Exact match $genSeqID:$genStart-" . ($genStart+$matchLen-1) . "_" . $orient . "\n";
          $matches[$seqRecsIdx] = [] if ( ! defined $matches[$seqRecsIdx] );
          push @{$matches[$seqRecsIdx]},[$genSeqID,$genStart,$genStart+$matchLen-1,$orient];
        }
      }
    }
  }
  close IN;

  # Extract sequences from the reference 2bit file
  #open IN,"$pathToUCSCTools/twoBitToFa -bed=$tmpFilename $twoBitRef stdout|" or 
  #   croak "validateSequences: Error could not open $twoBitRef with $pathToUCSCTools/twoBitToFa!\n";
  #my $tmp_seq;
  #my $seq = "";
  #my @refSeqs = ();
  #$idx = "";
  #while (<IN>){
  #  if ( />seq-(\d+)/ ) {
  #    my $tIdx = $1;
  #    $refSeqs[$idx] = uc($seq) if ( $seq ne "" );
  #    $seq = "";
  #    $idx = $tIdx;
  #    next;
  #  }
  #  s/[\n\r\s]+//g;
  #  $seq .= $_;
  #}
  #$refSeqs[$idx] = uc($seq) if ( $seq ne "" );
  #close IN;
  #unlink($tmpFilename);
  #
  ## Validate sequences/ranges against the reference set
  #my @validationErrors = ();
  #for ( my $i = 0; $i <= $#seqRanges; $i++ ){
  #  # Unpack datastructures
  #  my $idRec = $seqRanges[$i];  
  #  my $start = $idRec->[1];
  #  my $end = $idRec->[2];
  #  my $orient = $idRec->[3];
  #  # The length of the pads acheived on this particular 
  #  # sequence. This may be less than the padding specified
  #  # if the at the begining or end of a contig.
  #  my $leftPad = $leftPads[$i];
  #  my $rightPad = $rightPads[$i];
  #  # The original sequence and un-parsed identifier
  #  my $origID = $seqRecs->[$i]->[0];
  #  my $origSeq = $seqRecs->[$i]->[1];
  #  # and its reverse complement
  #  my $revOrigSeq = reverse($origSeq);
  #  $revOrigSeq =~ tr/ACGTRYWSKMNXBDHV/TGCAYRSWMKNXVHDB/;
  #
  #  # The reference sequence pulled from the twobit file
  #  my $refSeq = $refSeqs[$i];
  #
  #  # Identify matches by shifting over all possible
  #  # origSeq-length substrings of the padded reference
  #  # sequence.
  #  my $lenSeq = length($origSeq);
  #  my $lenDiff = length($refSeq) - $lenSeq;
  #  my @possibleMatches = ();
  #  my $startOffset = '';
  #  my $endOffset = '';
  #  my $orientStr = '';
  #  my $validMatch = 0;
  #  for ( my $j = 0; $j <= $lenDiff; $j++ ) {
  #    # Convert this $j index into offsets from the
  #    # original start and end sequence indexes.
  #    my $sOffset = $j - $leftPad;
  #    my $eOffset = $sOffset + ($lenSeq-($end-$start+1));
  #
  #    # Compare subsequences in both orientations
  #    my $matchDir = '';
  #    if ( $orient eq "+" ) {
  #      # Assume '+' first
  #      if ( $origSeq eq substr($refSeq,$j,$lenSeq) ) {
  #        $matchDir = '+'; # Match to positive strand
  #      }elsif ( $revOrigSeq eq substr($refSeq,$j,$lenSeq) ) {
  #        $matchDir = '-'; # Match to reverse strand
  #      }
  #    }else {
  #      # Assume '-' first
  #      if ( $revOrigSeq eq substr($refSeq,$j,$lenSeq) ) {
  #        $matchDir = '-'; # Match to reverse strand
  #      }elsif ( $origSeq eq substr($refSeq,$j,$lenSeq) ) {
  #        $matchDir = '-'; # Match to positive strand
  #      }
  #    }
  #
  #    # If we have a match
  #    if ( $matchDir ne '' ) {
  #      # If the match is where we were expecting to find it
  #      # (start offset == 0), process this as the preferred 
  #      # hit.
  #      if ( $sOffset == 0 ) {
  #        # Correct-ish.  It could be in the reverse, 
  #        # orientation, or at a different end offset.
  #        $startOffset = $sOffset;
  #        $endOffset = $eOffset;
  #        if ( $matchDir ne $orient ) {
  #          if ( $orient eq "+" ) {
  #            $orientStr = "-";
  #          }else {
  #            $orientStr = "+";
  #          }
  #        }else {
  #          # This flags that this was a perfect match
  #          $validMatch = 1 if ( $eOffset == 0 );
  #        }
  #      }else {
  #        # If there was a match and it's not where we expected
  #        # (start offset == 0), put this on a list of alternative
  #        # matches.
  #        push @possibleMatches, [$sOffset,$eOffset,$matchDir];
  #      }
  #    }
  #  }
  #  next if ( $validMatch );
  #
  #  if ( $startOffset ne '' ) {
  #    push @validationErrors, [$i, $origID, $startOffset, $endOffset, $orientStr, []];
  #  }elsif ( @possibleMatches ) {
  #    if ( scalar(@possibleMatches) == 1 ) {
  #      my $orientStr = '';
  #      if ( $possibleMatches[0]->[2] ne $orient ) {
  #        $orientStr = $possibleMatches[0]->[2];
  #      }
  #      push @validationErrors, [$i, $origID, $possibleMatches[0]->[0], $possibleMatches[0]->[1], $orientStr, []];
  #    }else { 
  #      push @validationErrors, [$i, $origID, "", "", "", \@possibleMatches];
  #    }
  #  }else {
  #    push @validationErrors, [$i, $origID, "", "", "", []];
  #  }
  #}
  #
  #return ( \@validationErrors );
}

 
1;
