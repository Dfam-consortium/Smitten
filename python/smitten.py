# -*- coding: utf-8 -*-
"""
smitten - Reference implementation of the Smitten sequence identifier format

        my ( $assemblyID, $sequenceID, $rangesRef ) =
             smitten::parseID( "hg38:chr1:100-200_+", 0 );

        my $nID = smitten::normalizeID( "hg38:chr1:100-200_-:10-20_+" );

        Library of functions to manipulate sequence identifiers in the Smitten
        format.  This includes legacy versions of the formats used by Arian
        Smit in his various tools.

        For example to normalize a sequence identifier in legacy, V1, or
        V2 Smitten format:

        legacy:
              chr1:1-200            forward strand range, seq identifier
              chr1:200-1            reverse strand range, seq identifier

        Smitten V1:
              chr1_100_200          positive strand range
              hg38:chr1_100_200     positive strand range with assembly ID
              hg38:chr1_100_200_R   negative strand range

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

        and the chained identifier:
             my ( $assemblyID, $sequenceID, $rangesRef ) =
                     smitten::parseID( "chr1:100-200_+:10-30_-" );

        would generate the following ranges structure:
                             [
                               [ 10, 30, '-']
                               [ 100, 200, '+' ],
                             ]

        Finally chained identifiers may also be normalized with the normalizeID() function:

           my $nID = smitten::normalizeID( "hg38:chr1:100-200_-:10-20_+" );

        which represents the positions 10 to 20 in the reversed sequence from
        100 to 200.  This routine would return "hg38:chr1:180-190_-".


    The Dfam project has settled on using the Black
    styleguide ( https://black.readthedocs.io ).

SEE ALSO: related_module.py
          Dfam: http://www.dfam.org

AUTHOR(S):
    Robert Hubley <rhubley@isbscience.org>

LICENSE:
    This code may be used in accordance with the Creative Commons
    Zero ("CC0") public domain dedication:
    https://creativecommons.org/publicdomain/zero/1.0/

DISCLAIMER:
  This software is provided ``AS IS'' and any express or implied
  warranties, including, but not limited to, the implied warranties of
  merchantability and fitness for a particular purpose, are disclaimed.
  In no event shall the authors or the Dfam consortium members be
  liable for any direct, indirect, incidental, special, exemplary, or
  consequential damages (including, but not limited to, procurement of
  substitute goods or services; loss of use, data, or profits; or
  business interruption) however caused and on any theory of liability,
  whether in contract, strict liability, or tort (including negligence
  or otherwise) arising in any way out of the use of this software, even
  if advised of the possibility of such damage.

"""

import re


def parseID( id_str, is_zerobased_halfopen=False ):
    """
    parseID() - blah

    Use: import smitten as sm
         ( assembly_id, sequence_id, ranges_ref ) =
                 sm.parseID( id_str, is_zerobased_halfopen )

    id_str               : Sequence identifier in V1 or V2 Smitten format or
                            legacy format (see below)

    is_zerobased_halfopen: For conversion of legacy identifiers, this treats
                            the coordinate ranges as zero-based, half-open
                            coordinates.

    Ranges are optional and may be in legacy, V1 or V2 formats.

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
    """
    orig_id_str = id_str
    id_prefix = id_str
    id_suffix = ""
    ranges = []

    # or... if re.search("[\s\n\r]",orig_id_str):
    if any([c in orig_id_str for c in {" ","\n","\r"}]):
        raise ValueError ( "identifier \'" + orig_id_str + "\' contains a space" +
                           "or a line termination character!\n" )

    version = None

    range_pattern = re.compile("(([:_])(\d+)([-_])(\d+)((_)([R\+\-]))?)")
    m_rng = re.search(range_pattern,id_str)
    if m_rng:
        id_prefix = id_str[0:m_rng.start()]
        id_str = id_str[m_rng.start():]
        for m_iter in re.finditer(range_pattern,id_str):
            (m_full,m_pre,m_start,m_sep1,m_end,m_suf,m_sep2,m_ori) = m_iter.groups()
            if ( m_pre == ":" and m_sep1 == "-" and m_sep2 == "_" ):
                version = 2
            elif ( m_pre == "_" and m_sep1 == "_" ):
                version = 1
            else:
                version = 0
            orient = "+"
            start = int(m_start)
            end = int(m_end)
            if start > end:
                if version > 0:
                    raise ValueError( "V1/V2 identifiers \'" + orig_id_str + "\' must have increasing range order!\n" )
                orient = "-"
                tmp = end
                end = start
                start = tmp
            elif m_ori == "":
                orient = "+"
            if m_ori == "R":
                orient = "-"
            if ( not is_zerobased_halfopen and ( start < 1 or end < 1 )):
                raise ValueError( "Sequence range less than 1 (" + orig_id_str + ").\n" )
            if ( len(ranges) and  (start > (ranges[-1][1] - ranges[-1][0] + 1) or
                 end > (ranges[-1][1] - ranges[-1][0] + 1) )):
                raise ValueError( "Sequence sub-range outside parent range ( " + orig_id_str + " ).\n" )
            if is_zerobased_halfopen:
                # Convert back into one-based, fully closed
                start = start + 1;
            ranges.append([start, end, orient])
            if m_iter.end() < (len(id_str) - 1):
                id_suffix = id_str[m_iter.end()+1:]
            else:
                id_suffix = ""

    if id_suffix != "":
        raise ValueError("Identifier \'" + orig_id_str + "\' has a non-standard suffix!\n")

    ids = id_prefix.split(":")
    assembly_id = None
    sequence_id = None
    if ( len(ids) == 2 ):
        assembly_id = ids[0]
        sequence_id = ids[1]
    elif len(ids) == 1:
        sequence_id = ids[0]
    if ( sequence_id is None or sequence_id == "" ):
        raise ValueError("Identifier \'" + orig_id_str + "\' does not have a sequence identifier!\n")

    return( assembly_id, sequence_id, ranges, version )


#
#    def normalizeID(self):
#    """
#    normalizeID()
#
#  Use:  my $normalizedIDStr = smitten::normalizeID( $idStr );
#
#      $idStr : Sequence identifier in legacy, V1 or V2 Smitten format.
#
#  Returns the normalized smitten identifier consisting of a single
#  sequence range given an identifier that has more than one level of 
#  sequence ranges.  For example, given the input identifier:
#
#          hg38:chr1:100-200_-:10-20_+
#
#  representing the positions 10 to 20 in the reversed sequence from
#  100 to 200.  This routine would return:
#
#          hg38:chr1:180-190_-
#
#    """
#    my $idStr = shift;
#
#    my ( $assemblyID, $sequenceID, $rangesRef ) = parseID( $idStr );
#
#    # Nothing to normalize
#    return $idStr if ( @{$rangesRef} == 0 );
#
#    my $startIdx = $rangesRef->[$#{$rangesRef}]->[0];
#    my $endIdx = $rangesRef->[$#{$rangesRef}]->[1];
#    my $currOrient = $rangesRef->[$#{$rangesRef}]->[2];
#    for ( my $i = $#{$rangesRef}-1; $i >= 0; $i-- ) {
#      if ( $rangesRef->[$i]->[2] eq "-" ) {
#        $startIdx = $rangesRef->[$i]->[1]-$startIdx+1;
#        $endIdx = $rangesRef->[$i]->[1]-$endIdx+1;
#      }else {
#        $startIdx = $rangesRef->[$i]->[0]+$startIdx-1;
#        $endIdx = $rangesRef->[$i]->[0]+$endIdx-1;
#      }
#      if ( $rangesRef->[$i]->[2] eq '-' && $currOrient eq '-' ) {
#        $currOrient = '+';
#      }elsif ( $rangesRef->[$i]->[2] eq '-' || $currOrient eq '+' ) {
#      $currOrient = '-';
#      }
#    }
#
#    my $retStr = "";
#    $retStr .= $assemblyID . ":" if ( $assemblyID ne "" );
#    $retStr .= $sequenceID . ":";
#    if ( $startIdx < $endIdx ) {
#      $retStr .= $startIdx . "-" . $endIdx;
#    }else {
#      $retStr .= $endIdx . "-" . $startIdx;
#    }
#    $retStr .= "_$currOrient";
#
#    return $retStr;
#
###    }
