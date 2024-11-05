# Overview

This library provides reference implementations for parsing and normalizing 
sequence identifiers in the Smitten format. The format allows for the encoding
of DNA sequence ranges and strand orientation to existing sequence identifiers.
These ranges may be recursively defined by adding further subranges.  This 
format has been in-use in many sequence analysis tools/scripts developed by 
Arian Smit to process Transposable Element sequences.  The format has evolved 
over the years and now encompases two legacy version (V0/V1), and the current 
supported identifier format (V2).  

At this time there is a robust Rust and Perl reference implementation in the
rust and perl sudirectories.  A Python reference implementation is in-progress.


# Format Specification

There are three versions of the Smitten format that are supported by this
library.  Each format is described below first with a formal specification
in Augmented Backus-Naur Form (ABNF), and then with examples.

## V0 Format
```abnf
ID = sequence_identifier *[ “_” start_position “_” end_position [“_R”]]
special_chars = ":" / “/” / “?” / “#” / “[“ / “]” / “@” / “!” / “$” / “&” / “’” / “(“ / “)”
                “*” / “+” / “,” / “;” / “=“ / “~” / “|” / “^” / “”” / “>” / “<“ / “.” / “%”
                "-" / "_"
sequence_identifier = 1*(ALPHA / DIGIT / special_chars)
start_position = 1*(DIGIT) ; 1-based sequence position
end_position = 1*(DIGIT) ; 1-based sequence position, fully closed coordinates
```

```plaintext
Examples:                          Meaning:
Chr1                               seq_id="Chr1", whole chromosome
1_10_30                            seq_id="1", from 10-30, forward strand     
seq1_50_100                        seq_id="seq1", from 50-100, forward strand
seq1_1_10_30                       seq_id="seq1_1", from 10-30, forward strand
seq1:2_10_30                       seq_id="seq1:2", from 10-30, forward strand
seq1_1_100_10_30_R                 seq_id="seq1", from 10-30 reverse strand of 1-100 forward strand
chr1_11023_38232_R_100_200         seq_id="chr1", from 100-200 forward strand of 11023-38232 reverse strand   
seq1_exon2_100_200_R               seq_id="seq1_exon2", from 100-200 reverse strand
```

Pitfalls:
This format relies on the suffix sequence always containing a meaningful set of range/orientation tokens
separated by a non-reserved underscore ("_") character.  Identifiers ending in a "_#_#" pattern are 
indistinguishable from a subrange specification.

## V1 Format
```abnf
ID = sequence_identifier *[ “:” (forward_orient-range / reverse_orient-range) ]
special_chars = ":" / “/” / “?” / “#” / “[“ / “]” / “@” / “!” / “$” / “&” / “’” / “(“ / “)”
                “*” / “+” / “,” / “;” / “=“ / “~” / “|” / “^” / “”” / “>” / “<“ / “.” / “%”
                "-" / "_"
reverse_orient_range = upper_bound “-” lower_bound
forward_orient_range = lower_bound “-” upper_bound
lower_bound = 1*(DIGIT) ; 1-based sequence position
upper_bound = 1*(DIGIT) ; 1-based sequence position, fully closed coordinates
```

Examples:                          Meaning:
```
Examples:
Chr1                               seq_id="Chr1", whole chromosome
Seq1:10-30                         seq_id="Seq1", from 10-30, forward strand
Seq1:30-10                         seq_id="Seq1", from 10-30, reverse strand
Seq1:100-200:10-30                 seq_id="Seq1", from 10-30 forward strand of 100-200 forward strand
Seq1:100-200:30-10                 seq_id="Seq1", from 10-30 reverse strand of 100-200 forward strand
```

Pitfalls:
The use of range ordering to denote orientation is problematic in a 1-based coordinate system. The
strand of a single base position range is ambiguous.


## V2 Format
```abnf
ID = [assembly_identifier “:”] sequence_identifier *[“:” start_position “-” end_position “_” orient]
special_chars = “/” / “?” / “#” / “[“ / “]” / “@” / “!” / “$” / “&” / “’” / “(“ / “)”
                “*” / “+” / “,” / “;” / “=“ / “~” / “|” / “^” / “”” / “>” / “<“ / “.” / “%”
                "-" / "_"
orient = “+” / “-”
assembly_identifier = *(ALPHA / DIGIT / special_chars )
sequence_identifier = 1*(ALPHA / DIGIT / special_chars )
start_position = 1*(DIGIT) ; 1-based sequence position
end_position = 1*(DIGIT) ; 1-based sequence position, fully closed coordinates
```

```plaintext
Examples:                          Meaning:
Chr1
seq_1:30-40_+
Seq1;contig4:100-103_-
chr1:11023-38232_-:100-200_+
hg38:chr1
GCA30283222:contig_1:100-200_+
```


-Robert Hubley, 2022-2024
