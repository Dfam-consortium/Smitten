#!/usr/bin/perl
use FindBin;
use lib $FindBin::RealBin;
use smitten;
use Data::Dumper;

my $VERBOSE = 0;
if ( $ARGV[0] ) {
  $VERBOSE = 1;
}


my @convertRecs = (
##
## V0 Examples
## 
## Typical examples
# test_string              coord_type parse_fmt exp_outcome exp_version exp_conv
["chr1_100_200",              "obfc",   "any",    "pass",       0,       "chr1:100-200_+"],
["seq_100_200_R",             "obfc",   "any",    "pass",       0,       "seq:100-200_-"],
["(ACC)n#Simple_1_10",        "obfc",   "any",    "pass",       0,       "(ACC)n#Simple:1-10_+"],
["chr_1_100_200_10_20_R",     "obfc",   "any",    "pass",       0,       "chr_1:100-200_+:10-20_-"],
["AMM1015_1000_2000_R_1_10_R","obfc",   "any",    "pass",       0,       "AMM1015:1000-2000_-:1-10_-"],
## Atypical but will be recognized as V0 and converted to V2
["hg38:chr1_100_200",         "obfc",   "any",    "pass",       0,       "hg38:chr1:100-200_+"],
["hg38:chr1_100_200_R",       "obfc",   "any",    "pass",       0,       "hg38:chr1:100-200_-"],
["1_2000_3000_400_500_60_70_8_9", "obfc","any",   "pass",       0,       "1:2000-3000_+:400-500_+:60-70_+:8-9_+"],
["chr;1_100_200",             "obfc",   "any",    "pass",       0,       "chr;1:100-200_+"],
## Incorrect V0 cases
# Sequence subrange outside range
["chr_1_100_200_150_200_R",   "obfc",   "any",    "fail",      "",       ""],
["HSPA2_0_0",                 "obfc",   "any",    "fail",      "",       ""],
# Mixture of V0/V1 concepts (V0 with reverse index ordering)
["chr2_200_100",              "obfc",   "any",    "fail",      "",       ""],
["NT_004321_0",               "obfc",   "any",    "fail",      "",       ""],
["chr13:51174549-51174548_R", "obfc",   "any",    "fail",      "",       ""],
# Missing identifier
["_100_200",                "obfc", "any", "fail", "",""],
## V0 cases that cannot be converted to V2
# The assembly/sequence sequence identifier would contain more than one ":" character
["AMM::1002:Seq1:100_200",    "obfc",   "any",    "fail",      "",       ""],
# The ":" is a valid sequence identifier in V0, but not in V2
[":_100_200",                  "obfc",   "any",    "fail",      "",       ""],

##
## V1 Examples
## 
## Typical examples (found in old Dfam seed alignment files)
["chr1:1-200",                 "obfc", "any", "pass", 1, "chr1:1-200_+"],
["chr1:200-1",                 "obfc", "any", "pass", 1, "chr1:1-200_-"],
["GCA1:chr2:540133-28232",     "obfc", "any", "pass", 1, "GCA1:chr2:28232-540133_-"],
# Multi-range format -- less common
["chr1:200-1:10-5",            "obfc", "any", "pass", 1, "chr1:1-200_-:5-10_-"],
["(ACC)n#Simple:400-500:20-1", "obfc", "any", "pass", 1, "(ACC)n#Simple:400-500_+:1-20_-"],
## Atypical examples
["seq1@*@)(_:1-200",              "obfc", "any", "pass",   1, "seq1@*@)(_:1-200_+"],
["1:2000-3000:400-500:60-70:8-9", "obfc", "any", "pass",   1, "1:2000-3000_+:400-500_+:60-70_+:8-9_+"],
["PB2_4_8_8:2709-6535",           "obfc", "any", "pass",   1, "PB2_4_8_8:2709-6535_+"],
["chromosome___:0-10",            "zbho", "any", "pass",   1, "chromosome___:1-10_+"],
## Incorrect V1 cases
# Sequence subrange outside range
["chr_1:100-200:200-150",    "obfc",   "any",    "fail",      "",       ""],
# zero indices (with obfc)
["ASeq:0-0",                 "obfc",   "any",    "fail",      "",       ""],
# Mixture of V2/V1/V0 concepts
["chr1:200-1:10-5_+",       "obfc",   "any",    "fail",       "",       ""],
["chr1:200-1_+",            "obfc",   "any",    "fail",       "",       ""],
["hg38:chr1:200-100_-",     "obfc",   "any",    "fail",       "",       ""],
## V1 cases that cannot be converted to V2
# The assembly/sequence sequence identifier would contain more than one ":" character
["AMM::1002:Seq1:100-200",    "obfc",   "any",    "fail",      "",       ""],
# The ":" is a valid sequence identifier in V1, but not in V2
["::100-200",                "obfc", "any", "fail", "",""],


##
## V2 Examples
##
## Typical examples
["chr1:100-200_+",          "obfc",   "any",    "pass",       2,       "chr1:100-200_+"],
["chr1_1_5_3:1-3_+",        "obfc", "any", "pass", 2 , "chr1_1_5_3:1-3_+"],
["chr1:100-200_+",          "obfc", "any","pass",2,"chr1:100-200_+"],     
["hg38:chr1:100-200_+",     "obfc", "any","pass", 2,"hg38:chr1:100-200_+"], 
["hg38:chr1:100-200_-",     "obfc", "any","pass",2,"hg38:chr1:100-200_-"],
["hg38:chr1:100-200_+:10-50_-:1-5_+", "obfc", "any","pass",2,"hg38:chr1:100-200_+:10-50_-:1-5_+"],  
["hg38:chr1:200-200_-",     "obfc", "any","pass",2,"hg38:chr1:200-200_-"],
["hg38:chr1:33438223-33439283_+", "zbho", "any","pass",2,"hg38:chr1:33438224-33439283_+"],
["JANCRE010000006.1:5563658-5564462_-:1-458_+", "obfc", "any", "pass", 2,"JANCRE010000006.1:5563658-5564462_-:1-458_+"],
["JANCRE010000006.1_5563658_5564462_R:1-458_+", "obfc", "any", "pass", 2,"JANCRE010000006.1_5563658_5564462_R:1-458_+"],
["hg_38:chr+1:10-40_+",     "obfc", "any", "pass", 2,"hg_38:chr+1:10-40_+"],
["hg-38:chr_1:10-40_+",     "obfc", "any", "pass", 2,"hg-38:chr_1:10-40_+"],

##
## Unrecognizable Examples
##
# Unallowed characters 
["seq1_ 1_2",               "obfc", "any", "fail", "",""],
["seq1:1-2\n_+",            "obfc", "any", "fail", "",""],
# Sequences without a range
["chr1",                    "obfc", "any", "pass", "","chr1"],
["hg38:chr1",               "obfc", "any", "pass", "","hg38:chr1"],
["100_200",                 "obfc", "any", "pass", "","100_200"],
["100-200",                 "obfc", "any", "pass", "","100-200"],
["100:200",                 "obfc", "any", "pass", "","100:200"],
# Failures based on incorrect ":"s
["100:200:",                 "obfc", "any", "fail", "",""],
["100:200:seq:",                 "obfc", "any", "fail", "",""],

);


if ( $VERBOSE ) {
  print "ConvertID Tests\n";
  print "---------------\n";
}
foreach my $testRec ( @convertRecs ) {
  my $idStr = $testRec->[0];
  my $coordType = $testRec->[1];
  my $parseFmt = $testRec->[2];
  my $expOutcome = $testRec->[3];
  my $expVersion = $testRec->[4];
  my $expConv = $testRec->[5];
  my $failDesc = "";

  my $isZeroBasedHalfOpen = 0;
  $isZeroBasedHalfOpen = 1 if ( $coordType eq "zbho" );
 
  print "Testing convertID($idStr, $isZeroBasedHalfOpen):\n" if ( $VERBOSE );
  my $newID, $vers;
  #eval {( $newID, $vers ) = convertID( id => $idStr, zbho => $isZeroBasedHalfOpen, format => $parseFmt );};
  # format parameter is deprecated
  eval {( $newID, $vers ) = convertID( id => $idStr, zbho => $isZeroBasedHalfOpen );};
  if ( $@ ) {
    if ( $expOutcome ne "pass" ){
      print "  Ok.\n" if ( $VERBOSE );
      next;
    }
    $failDesc .= "$idStr: unexpected conversion failure $@\n";
  }elsif ( $expOutcome eq "fail" ) {
    $failDesc .= "$idStr: unexpected parsing success\n";
  }
  $failDesc .= "$idStr: failed conversion: received $newID expected $expConv\n" 
    if ( $newID ne $expConv );
  $failDesc .= "$idStr: failed version identification: received $vers expected $expVersion\n" 
    if ( $vers != $expVersion );

  if ( $failDesc ne "" ) {
    print $failDesc;
  }elsif ( $VERBOSE ) {
    print "  Ok.\n";
  }
}



# Test Encoding:
#   Test String: The complete ID string to be parsed in the test
#   Expected Outcome: The expected parsability of the test string either
#                     "pass" or "fail".
#   Expected Asssembly: The expected returned assembly string (if present)
#                       otherwise "".
#   Expected Sequence ID: The expected returned sequence identifier.
#   Epected Ranges: A reference to a list of list references containing the
#                   expected start, end and orientation.
#   Expected Version: The expected identifier format version, where
#                     0 = V0, 1 = V1, and 2 = V2.
#   Expected Normalisation: The expected normalized range value.
#                         
my @parseRecs = ( 
["chr1:100-200_+", "pass","", "chr1", [[100, 200, '+']],2,"chr1:100-200_+"],     
["hg38:chr1:100-200_+", "pass","hg38", "chr1", [[100, 200, '+']],2,"hg38:chr1:100-200_+"], 
["hg38:chr1:100-200_-", "pass","hg38", "chr1", [[100, 200, '-']],2,"hg38:chr1:100-200_-"],
["hg38:chr1:100-200_+:10-50_-:1-5_+", "pass","hg38", "chr1", [[100, 200, '+'],[10,50,'-'],[1,5,'+']],2,"hg38:chr1:145-149_-"],  
["hg38:chr1:200-200_-", "pass","hg38", "chr1", [[200, 200, '-']],2,"hg38:chr1:200-200_-"],
["hg38:chr1:200-100_-", "fail","", "", [[]],"",""],
["hg38:chr1:33438223-33439283_+", "pass","hg38", "chr1", [[33438223,33439283,'+']],2,"hg38:chr1:33438223-33439283_+"],
["JANCRE010000006.1:5563658-5564462_-:1-458_+", "pass", "", "JANCRE010000006.1", [[5563658,5564462,'-'],[1,458,'+']],2,"JANCRE010000006.1:5564005-5564462_-"],
["hg_38:chr+1:10-40_+", "pass", "hg_38", "chr+1", [[10,40,'+']],2,"hg_38:chr+1:10-40_+"],
["hg-38:chr_1:10-40_+", "pass", "hg-38", "chr_1", [[10,40,'+']],2,"hg-38:chr_1:10-40_+"],
["chr1_1_5_3:1-3_+", "pass", "", "chr1_1_5_3", [[1,3,'+']],2,"chr1_1_5_3:1-3_+"],
["scaffold_164_5", "pass", "", "scaffold_164_5", [], 2, "scaffold_164_5"],
["hg38:mm10:chr1:100-200_-", "fail", "", "", [], "","" ],
#   V2   seqid=scaffold_164_5
);

if ( $VERBOSE ) {
  print "parseID and normalizeID Tests\n";
  print "-----------------------------\n";
}
foreach my $testRec ( @parseRecs ) {
  my $idStr = $testRec->[0];
  my $expOutcome = $testRec->[1];
  my $expAssembly = $testRec->[2];
  my $expSeqID = $testRec->[3];
  my $expRanges = $testRec->[4];
  my $expVersion = $testRec->[5];
  my $expNormalisation = $testRec->[6];
  my $failDesc = "";

  print "Testing parseID($idStr):\n" if ( $VERBOSE );
  my $a, $s, $r;
  eval {( $a, $s, $r ) = parseID( $idStr );};
  if ( $@ ) {
    if ( $expOutcome ne "pass" ) {
      print "  Ok.\n" if ( $VERBOSE );
      next;
    }
    $failDesc .= "$idStr: unexpected parsing failure $@\n";
  }elsif ( $expOutcome eq "fail" ) {
    $failDesc .= "$idStr: unexpected parsing success\n";
  }
  $failDesc .= "$idStr: failed assembly parse: received $a expected $expAssembly\n" 
    if ( $a ne $expAssembly );
  $failDesc .= "$idStr: failed sequence ID parse: received $s expected $expSeqID\n" 
    if ( $s ne $expSeqID );
  if ( @{$expRanges} != @{$r} ) {
    $failDesc .= "$idStr: failed, range count mismatch! Expected:\n" . Dumper($expRanges) . "\nReceived:\n" . Dumper($r) . "\n"
  }else {
    for ( my $i = 0; $i <= $#{$expRanges}; $i++ ) {
      if ( @{$expRanges->[$i]} != @{$r->[$i]} || 
           $expRanges->[$i]->[0] != $r->[$i]->[0] ||
           $expRanges->[$i]->[1] != $r->[$i]->[1] ||
           $expRanges->[$i]->[2] != $r->[$i]->[2] ) {
        $failDesc .= "$idStr: failed, range mismatch! Expected:\n" . Dumper($expRanges) . "\nReceived:\n" . Dumper($r) . "\n";
        last;
      }
    }
  }
  if ( $VERBOSE && $failDesc eq "" ) {
    print "  Ok.\n";
  }

  print "Testing normalizeID($idStr):\n" if ( $VERBOSE );
  my $n = normalizeID($idStr);
  #print "$idStr:  Normalized = $n\n";
  $failDesc .= "$idStr: failed normalisation: received $n expected $expNormalisation\n" 
    if ( $n ne $expNormalisation );

  if ( $failDesc ne "" ) {
    print $failDesc;
  }elsif ( $VERBOSE ) {
    print "  Ok.\n";
  }

}

#foreach my $idStr ( "hg38:chr1:100-200[+]", "hg38:chr1:100-200[-]:10-20[+]", "hg38:chr1_100_200", "hg38:chr1_100_200_R", "chr1_100_200_R_10_20_R", "chr1_100_200_R_10_20" ) {
#  my ( $a, $s, $r ) = parseID( $idStr );
#  print "$idStr : assembly = $a, sequence = $s, ranges = " . Dumper($r) . "\n";
#  print "Normalized = " . normalizeID($idStr) . "\n";
#}

