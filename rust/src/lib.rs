/*!
This library provides routines for parsing and normalizing sequence identifiers 
in the Smitten format. The format allows for the encoding of DNA sequence ranges 
and strand orientation to existing sequence identifiers.  These ranges may be
recursively defined by adding further subranges.  This format has been in-use
in many sequence analysis tools/scripts developed by Arian Smit to process 
Transposable Element sequences.  The format has evolved over the years and now
encompases two legacy version (V0/V1), and the current supported identifier 
format (V2).  

Here is a quick overview of how the library is used:

 ```rust
 use smitten::{Identifier, IDVersion};

 // Parse a current Smitten identifier
 let parsed_id = Identifier::from_v2("hg38:chr1:100-200_+").unwrap();
 assert_eq!(parsed_id.assembly_id, Some("hg38".to_string()));
 assert_eq!(parsed_id.sequence_id, "chr1".to_string());
 assert_eq!(parsed_id.ranges.len(), 1);
 assert_eq!(parsed_id.ranges[0].start, 100);
 assert_eq!(parsed_id.ranges[0].end, 200);
 assert_eq!(parsed_id.ranges[0].orientation, '+');

 // Parse an unknown Smitten identifier and convert to V2
 let (parsed_id, inferred_version) = Identifier::from_unknown_format("chr1_100_200", false).unwrap();
 assert_eq!(parsed_id.to_string(), "chr1:100-200_+");

 // Parse a multi-range Smitten identifier and normalize
 let parsed_id = Identifier::from_v2("hg38:chr1:100-200_+:10-50_-:1-5_+").unwrap();
 assert_eq!(parsed_id.ranges.len(), 3);
 let normalized_id = parsed_id.normalize().unwrap();
 assert_eq!(normalized_id.to_string(), "hg38:chr1:145-149_-");
 ```

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
Examples:
Chr1                               seq_id="Chr1", whole chromosome                                
1_10_30                            seq_id="1", from 10-30, forward strand
seq1_50_100                        seq_id="seq1", from 50-100, forward strand
seq1_1_10_30                       seq_id="seq1_1", from 10-30, forward strand
seq1:2_10_30                       seq_id="seq1:2", from 10-30, forward strand
seq1_1_100_10_30_R                 seq_id="seq1", from 10-30 reverse strand of 1-100 forward strand
chr1_11023_38232_R_100_200         seq_id="chr1", from 100-200 forward strand of 11023-38232 reverse strand
seq1_exon2_100_200_R               seq_id="seq1_exon2", from 100-200 reverse strand  

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

Examples:
Chr1                               seq_id="Chr1", whole chromosome
Seq1:10-30                         seq_id="Seq1", from 10-30, forward strand
Seq1:30-10                         seq_id="Seq1", from 10-30, reverse strand
Seq1:100-200:10-30                 seq_id="Seq1", from 10-30 forward strand of 100-200 forward strand
Seq1:100-200:30-10                 seq_id="Seq1", from 10-30 reverse strand of 100-200 forward strand

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
*/

use regex::Regex;

#[derive(Debug, PartialEq, Clone)]
pub enum IDVersion {
    Undefined,
    V0,
    V1,
    V2,
}

#[derive(Debug, Clone)]
pub struct Range {
    pub start: usize,
    pub end: usize,
    pub orientation: char, // '+' or '-'
}

#[derive(Debug, Clone)]
pub struct Identifier {
    pub assembly_id: Option<String>,
    pub sequence_id: String,
    pub ranges: Vec<Range>,
    pub inferred_version: IDVersion,
}

// Define the API
impl Identifier {
    /// Creates an `Identifier` from an identifier of unknown format (V0 or V1), converting it to V2.
    pub fn from_unknown_format(id: &str, zbho: bool) -> Result<(Self, IDVersion), String> {
        // Attempt to convert to V2 format
        let (v2_id, inferred_version) = Identifier::convert_id(id, zbho)?;

        // Parse the V2 identifier string into an `Identifier` struct
        let identifier = Identifier::parse_id(&v2_id)?;

        // Return both the parsed `Identifier` and the inferred version
        Ok((identifier, inferred_version))
    }

    // By providing these public APIs, we can provide more focused 
    // converters in the future for improved error handling.
    /// Other constructors for specific versions remain the same
    pub fn from_v0(id: &str) -> Result<Self, String> {
        let (v2_id, _) = Identifier::convert_id(id, false)?;
        Identifier::parse_id(&v2_id)
    }

    pub fn from_v1(id: &str) -> Result<Self, String> {
        let (v2_id, _) = Identifier::convert_id(id, false)?;
        Identifier::parse_id(&v2_id)
    }

    pub fn from_v2(id: &str) -> Result<Self, String> {
        Identifier::parse_id(id)
    }

    pub fn normalize(&self) -> Result<Self, String> {
        let normalized_id_str = self.normalize_id()?;
        Identifier::parse_id(&normalized_id_str)
    }
}


impl std::fmt::Display for IDVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let version_str = match self {
            IDVersion::Undefined => "Undefined",
            IDVersion::V0 => "V0",
            IDVersion::V1 => "V1",
            IDVersion::V2 => "V2",
        };
    write!(f, "{}", version_str)
    }
}

impl std::fmt::Display for Identifier {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let mut v2_id = String::new();

        if let Some(assembly) = &self.assembly_id {
            v2_id.push_str(&format!("{}:", assembly));
        }

        v2_id.push_str(&self.sequence_id);

        if !self.ranges.is_empty() {
            for range in &self.ranges {
                v2_id.push_str(&format!(":{}-{}_{}", range.start, range.end, range.orientation));
            }
        }

        write!(f, "{}", v2_id)
    }
}

impl Identifier {
    /// Converts a Smitten format identifier in V0, V1, or V2 to V2 format.
    ///
    /// # Arguments
    ///
    /// * `id` - Sequence identifier in V0, V1, or V2 format
    /// * `zbho` - Boolean flag; if true, treats coordinate ranges as zero-based half-open (ZBHO)
    ///
    /// # Returns
    ///
    /// Returns a tuple `(String, IDVersion)` with the V2 equivalent identifier and inferred version.
    ///
    fn convert_id(id: &str, zbho: bool) -> Result<(String, IDVersion), String> {
        if id.contains(|c: char| c.is_whitespace() || c == '\n' || c == '\r') {
            return Err(format!(
                "convertID: Identifier '{}' contains a space or a line termination character!",
                id
            ));
        }

        let re = Regex::new(r"(.*)(([:_])(\d+)([-_])(\d+)((_)([R+\-]))?)$").unwrap();
        let mut inferred_fmt = None;
        let mut sequence_id = id.to_string();
        let mut ranges = Vec::new();

        // Iterate over matches and parse ranges
        while let Some(captures) = re.captures(&sequence_id) {
            let start = captures[4].parse::<usize>().unwrap();
            let end = captures[6].parse::<usize>().unwrap();
            let orientation = captures.get(9).map_or('+', |m| {
                if m.as_str() == "R" { '-' } else { m.as_str().chars().next().unwrap() }
            });

            let start = if zbho { start + 1 } else { start };

            // Infer format based on separators if not already inferred
            let range_fmt = if captures.get(3).map(|s| s.as_str()) == Some(":")
                && captures.get(5).map(|s| s.as_str()) == Some("-")
                && captures.get(7).is_some()
            {
                IDVersion::V2
            } else if captures.get(3).map(|s| s.as_str()) == Some(":")
                && captures.get(5).map(|s| s.as_str()) == Some("-")
            {
                IDVersion::V1
            } else {
                IDVersion::V0
            };

            if let Some(ref inferred) = inferred_fmt {
                if *inferred != range_fmt {
                    break;
                }
            } else {
                inferred_fmt = Some(range_fmt.clone());
            }

            let (ordered_start, ordered_end, final_orientation) = match range_fmt {
                IDVersion::V0 | IDVersion::V2 if start > end => {
                    return Err(format!(
                        "convertID: V{} identifier '{}' must have increasing range order!",
                        range_fmt, sequence_id
                    ));
                }
                IDVersion::V1 if start > end => (end, start, '-'),
                IDVersion::V1 => (start, end, '+'),
                _ => (start, end, orientation),
            };

            ranges.push(Range {
                start: ordered_start,
                end: ordered_end,
                orientation: final_orientation,
            });

            sequence_id = captures[1].to_string();
        }

        if ranges.is_empty() {
            inferred_fmt = Some(IDVersion::Undefined);
        }

        let ids: Vec<&str> = sequence_id.split(':').collect();
        let (assembly_id, sequence_id) = match (sequence_id.matches(':').count(), ids.len()) {
            (1, 2) if !ids[0].is_empty() && !ids[1].is_empty() => (Some(ids[0]), ids[1]),
            (0, 1) if !sequence_id.is_empty() => (None, sequence_id.as_str()),
            _ => {
                return Err(format!(
                    "convertID: Identifier '{}' contains an invalid assembly+sequence structure.",
                    id
                ));
            }
        };

        let mut v2_id = if let Some(assembly) = assembly_id {
            format!("{}:{}", assembly, sequence_id)
        } else {
            sequence_id.to_string()
        };

        let mut current_parent_length = None;
        for range in ranges.iter().rev() {
            if range.start == 0 || range.end == 0 {
                return Err(format!(
                    "convertID: Invalid range {}-{} in a one-based fully-closed coordinate system.",
                    range.start, range.end
                ));
            }
            if let Some(parent_len) = current_parent_length {
                if range.start > parent_len || range.end > parent_len {
                    return Err(format!(
                        "convertID: Sequence sub-range {}-{} is outside the bounds of the parent range length {}.",
                        range.start, range.end, parent_len
                    ));
                }
            }
            current_parent_length = Some(range.end - range.start + 1);
            v2_id.push_str(&format!(":{}-{}_{}", range.start, range.end, range.orientation));
        }

        Ok((v2_id, inferred_fmt.unwrap_or(IDVersion::Undefined)))
    }

    /// Parses a Smitten format sequence identifier, returning the individual components.
    ///
    /// # Arguments
    ///
    /// * `id` - Sequence identifier in V0, V1, or V2 format
    ///
    /// # Returns
    ///
    /// Returns a `Identifier` struct containing `assembly_id`, `sequence_id`, `ranges`, and inferred `version`.
    ///
    fn parse_id(id: &str) -> Result<Self, String> {
        let re = Regex::new(r"(.*)(([:])(\d+)([-])(\d+)((_)([+\-]))?)$").unwrap();
        let mut assembly_id = None;
        let mut ranges = Vec::new();
        let mut id_str = id.to_string();

        // Remove ranges from the end of the ID string
        while let Some(captures) = re.captures(&id_str) {
            let start = captures[4].parse::<usize>().unwrap();
            let end = captures[6].parse::<usize>().unwrap();
            let orientation = captures[9].chars().next().unwrap();

            if start > end {
                return Err(format!(
                    "parseID: V2 identifiers '{}' must have increasing range order!",
                    id
                ));
            }

            ranges.push(Range { start, end, orientation });
            id_str = captures[1].to_string();
        }

        let ids: Vec<&str> = id_str.split(':').collect();
        let sequence_id;
        if ids.len() == 2 {
            assembly_id = Some(ids[0].to_string());
            sequence_id = ids[1].to_string();
        } else if ids.len() == 1 {
            sequence_id = ids[0].to_string();
        } else {
            return Err(format!(
                "parseID: Identifier '{}' contains an invalid number of ':' characters.",
                id
            ));
        }

        if sequence_id.is_empty() {
            return Err(format!(
                "parseID: Identifier '{}' does not have a sequence identifier!",
                id
            ));
        }

        Ok(Identifier {
            assembly_id,
            sequence_id,
            ranges: ranges.into_iter().rev().collect(),
            inferred_version: IDVersion::V2,
        })
    }

    /// Normalizes a chained sequence identifier in Smitten V2 format to a single normalized range.
    ///
    /// # Returns
    ///
    /// Returns a normalized sequence identifier as a `String`.
    ///
    fn normalize_id(&self) -> Result<String, String> {
        if self.ranges.is_empty() {
            let mut ret_str = String::new();
            if let Some(assembly) = &self.assembly_id {
                ret_str.push_str(&format!("{}:", assembly));
            }
            ret_str.push_str(&self.sequence_id);
            return Ok(ret_str);
        }

        let mut start_idx = self.ranges.last().unwrap().start;
        let mut end_idx = self.ranges.last().unwrap().end;
        let mut curr_orient = self.ranges.last().unwrap().orientation;

        for range in self.ranges.iter().rev().skip(1) {
            if range.orientation == '-' {
                start_idx = range.end - start_idx + 1;
                end_idx = range.end - end_idx + 1;
            } else {
                start_idx = range.start + start_idx - 1;
                end_idx = range.start + end_idx - 1;
            }

            curr_orient = if range.orientation == '-' && curr_orient == '-' {
                '+'
            } else if range.orientation == '-' || curr_orient == '+' {
                '-'
            } else {
                curr_orient
            };
        }

        let mut ret_str = String::new();
        if let Some(assembly) = &self.assembly_id {
            ret_str.push_str(&format!("{}:", assembly));
        }
        ret_str.push_str(&self.sequence_id);
        ret_str.push(':');

        if start_idx < end_idx {
            ret_str.push_str(&format!("{}-{}_{}", start_idx, end_idx, curr_orient));
        } else {
            ret_str.push_str(&format!("{}-{}_{}", end_idx, start_idx, curr_orient));
        }

        Ok(ret_str)
    }
}



#[cfg(test)]
mod tests {
    use super::*;

    fn run_convert_id_test(
        id: &str,
        coord_type: &str,
        exp_outcome: &str,
        exp_version: Option<IDVersion>,
        exp_v2_format: &str,
    ) {
        let zbho = coord_type == "zbho";
        let result = Identifier::from_unknown_format(id, zbho);

        match (result.clone(), exp_outcome, exp_version.clone()) {
            // Expected pass case
            (Ok((converted_id, version)), "pass", Some(expected_version)) => {
                assert_eq!(version, expected_version, "IDVersion mismatch for ID: {}", id);
                assert_eq!(converted_id.to_string(), exp_v2_format, "V2 format mismatch for ID: {}", id);
            }
            // Expected failure with None as version
            (Err(err), "fail", None) => {
                println!("Expected failure for ID: {}. Error message: {}", id, err);
            }
            // Expected failure with Some(Undefined) as version
            (Err(err), "fail", Some(IDVersion::Undefined)) => {
                println!("Expected specific failure for ID: {}. Error message: {}", id, err);
            }
            // Unexpected success when failure was expected
            (Ok((converted_id, version)), "fail", _) => {
                panic!(
                    "Unexpected success for ID: {}. Expected outcome: fail, but got success with V2 ID: '{}' and version: {:?}",
                    id, converted_id, version
                );
            }
            // Unexpected failure when pass was expected
            (Err(err), "pass", Some(expected_version)) => {
                panic!(
                    "Unexpected failure for ID: {}. Expected outcome: pass, version: {:?}. Error: {}",
                    id, expected_version, err
                );
            }
            // Catch-all for any unexpected outcome
            _ => panic!(
                "Unexpected outcome for ID: {}. Expected outcome: {}, version: {:?}, but got result: {:?}",
                id, exp_outcome, exp_version, result
            ),
        }
    }

    #[test]
    fn test_convert_id() {
        let test_cases = vec![
            // V0 Examples
            ("chr1_100_200", "obfc", "pass", Some(IDVersion::V0), "chr1:100-200_+"),
            ("seq_100_200_R", "obfc", "pass", Some(IDVersion::V0), "seq:100-200_-"),
            ("(ACC)n#Simple_1_10", "obfc", "pass", Some(IDVersion::V0), "(ACC)n#Simple:1-10_+"),
            ("chr_1_100_200_10_20_R", "obfc", "pass", Some(IDVersion::V0), "chr_1:100-200_+:10-20_-"),
            ("AMM1015_1000_2000_R_1_10_R", "obfc", "pass", Some(IDVersion::V0), "AMM1015:1000-2000_-:1-10_-"),
            ("hg38:chr1_100_200", "obfc", "pass", Some(IDVersion::V0), "hg38:chr1:100-200_+"),
            ("hg38:chr1_100_200_R", "obfc", "pass", Some(IDVersion::V0), "hg38:chr1:100-200_-"),
            ("1_2000_3000_400_500_60_70_8_9", "obfc", "pass", Some(IDVersion::V0), "1:2000-3000_+:400-500_+:60-70_+:8-9_+"),
            ("chr;1_100_200", "obfc", "pass", Some(IDVersion::V0), "chr;1:100-200_+"),
            ("chr_1_100_200_150_200_R", "obfc", "fail", Some(IDVersion::Undefined), ""),
            ("HSPA2_0_0", "obfc", "fail", Some(IDVersion::Undefined), ""),
            ("chr2_200_100", "obfc", "fail", Some(IDVersion::Undefined), ""),
            ("NT_004321_0", "obfc", "fail", Some(IDVersion::Undefined), ""),
            ("chr13:51174549-51174548_R", "obfc", "fail", Some(IDVersion::Undefined), ""),
            ("_100_200", "obfc", "fail", Some(IDVersion::Undefined), ""),
            ("AMM::1002:Seq1:100_200", "obfc", "fail", Some(IDVersion::Undefined), ""),
            (":_100_200", "obfc", "fail", Some(IDVersion::Undefined), ""),
            // V1 Examples
            ("chr1:1-200", "obfc", "pass", Some(IDVersion::V1), "chr1:1-200_+"),
            ("chr1:200-1", "obfc", "pass", Some(IDVersion::V1), "chr1:1-200_-"),
            ("GCA1:chr2:540133-28232", "obfc", "pass", Some(IDVersion::V1), "GCA1:chr2:28232-540133_-"),
            ("chr1:200-1:10-5", "obfc", "pass", Some(IDVersion::V1), "chr1:1-200_-:5-10_-"),
            ("(ACC)n#Simple:400-500:20-1", "obfc", "pass", Some(IDVersion::V1), "(ACC)n#Simple:400-500_+:1-20_-"),
            ("seq1@*@)(_:1-200", "obfc", "pass", Some(IDVersion::V1), "seq1@*@)(_:1-200_+"),
            ("1:2000-3000:400-500:60-70:8-9", "obfc", "pass", Some(IDVersion::V1), "1:2000-3000_+:400-500_+:60-70_+:8-9_+"),
            ("PB2_4_8_8:2709-6535", "obfc", "pass", Some(IDVersion::V1), "PB2_4_8_8:2709-6535_+"),
            ("chromosome___:0-10", "zbho", "pass", Some(IDVersion::V1), "chromosome___:1-10_+"),
            ("chr_1:100-200:200-150", "obfc", "fail", Some(IDVersion::Undefined), ""),
            ("ASeq:0-0", "obfc", "fail", Some(IDVersion::Undefined), ""),
            ("chr1:200-1:10-5_+", "obfc", "fail", Some(IDVersion::Undefined), ""),
            ("chr1:200-1_+", "obfc", "fail", Some(IDVersion::Undefined), ""),
            ("hg38:chr1:200-100_-", "obfc", "fail", Some(IDVersion::Undefined), ""),
            ("AMM::1002:Seq1:100-200", "obfc", "fail", Some(IDVersion::Undefined), ""),
            ("::100-200", "obfc", "fail", Some(IDVersion::Undefined), ""),
            // V2 Examples
            ("chr1:100-200_+", "obfc", "pass", Some(IDVersion::V2), "chr1:100-200_+"),
            ("chr1_1_5_3:1-3_+", "obfc", "pass", Some(IDVersion::V2), "chr1_1_5_3:1-3_+"),
            ("chr1:100-200_+", "obfc", "pass", Some(IDVersion::V2), "chr1:100-200_+"),
            ("hg38:chr1:100-200_+", "obfc", "pass", Some(IDVersion::V2), "hg38:chr1:100-200_+"),
            ("hg38:chr1:100-200_-", "obfc", "pass", Some(IDVersion::V2), "hg38:chr1:100-200_-"),
            ("hg38:chr1:100-200_+:10-50_-:1-5_+", "obfc", "pass", Some(IDVersion::V2), "hg38:chr1:100-200_+:10-50_-:1-5_+"),
            ("hg38:chr1:200-200_-", "obfc", "pass", Some(IDVersion::V2), "hg38:chr1:200-200_-"),
            ("hg38:chr1:33438223-33439283_+", "zbho", "pass", Some(IDVersion::V2), "hg38:chr1:33438224-33439283_+"),
            ("JANCRE010000006.1:5563658-5564462_-:1-458_+", "obfc", "pass", Some(IDVersion::V2), "JANCRE010000006.1:5563658-5564462_-:1-458_+"),
            ("JANCRE010000006.1_5563658_5564462_R:1-458_+", "obfc", "pass", Some(IDVersion::V2), "JANCRE010000006.1_5563658_5564462_R:1-458_+"),
            ("hg_38:chr+1:10-40_+", "obfc", "pass", Some(IDVersion::V2), "hg_38:chr+1:10-40_+"),
            ("hg-38:chr_1:10-40_+", "obfc", "pass", Some(IDVersion::V2), "hg-38:chr_1:10-40_+"),
            // Unrecognizable Examples
            ("seq1_ 1_2", "obfc", "fail", Some(IDVersion::Undefined), ""),
            ("seq1:1-2\n_+", "obfc", "fail", Some(IDVersion::Undefined), ""),
            ("chr1", "obfc", "pass", Some(IDVersion::Undefined), "chr1"),
            ("hg38:chr1", "obfc", "pass", Some(IDVersion::Undefined), "hg38:chr1"),
            ("100_200", "obfc", "pass", Some(IDVersion::Undefined), "100_200"),
            ("100-200", "obfc", "pass", Some(IDVersion::Undefined), "100-200"),
            ("100:200", "obfc", "pass", Some(IDVersion::Undefined), "100:200"),
            ("100:200:", "obfc", "fail", Some(IDVersion::Undefined), ""),
            ("100:200:seq:", "obfc", "fail", Some(IDVersion::Undefined), ""),
        ];

        for (id, coord_type, exp_outcome, exp_version, exp_v2_format) in test_cases {
            run_convert_id_test(id, coord_type, exp_outcome, exp_version, exp_v2_format);
        }
    }

fn run_parse_id_test(
    id: &str,
    expected_outcome: &str,
    expected_assembly: &str,
    expected_sequence_id: &str,
    expected_ranges: Vec<(usize, usize, char)>,
    expected_version: Option<IDVersion>,
    expected_normalization: &str,
) {
    let result = Identifier::from_v2(id);

    match (result.as_ref(), expected_outcome, expected_version.clone()) {
        // Expected pass case
        (Ok(parsed_id), "pass", Some(expected_version)) => {
            println!("Parsed ID: {:?}", parsed_id); // Debug output for parsed ID
            assert_eq!(
                parsed_id.assembly_id.as_deref().unwrap_or(""),
                expected_assembly,
                "Assembly mismatch for ID: {}", id
            );
            assert_eq!(parsed_id.sequence_id, expected_sequence_id, "Sequence ID mismatch for ID: {}", id);
            assert_eq!(parsed_id.inferred_version, expected_version, "IDVersion mismatch for ID: {}", id);

            let parsed_ranges: Vec<(usize, usize, char)> = parsed_id
                .ranges
                .iter()
                .map(|r| (r.start, r.end, r.orientation))
                .collect();
            assert_eq!(parsed_ranges, expected_ranges, "Range mismatch for ID: {}", id);

            // Normalization check if expected_normalization is provided
            if !expected_normalization.is_empty() {
                match parsed_id.normalize() {
                    Ok(normalized_id) => {
                        assert_eq!(normalized_id.to_string(), expected_normalization, "Normalization mismatch for ID: {}", id);
                    }
                    Err(err) => {
                        panic!("Normalization failed for ID: {}. Error: {}", id, err);
                    }
                }
            }



        }
        // Expected failure case
        (Err(err), "fail", None) => {
            println!("Expected failure for ID: {}. Error message: {}", id, err);
        }
        // Unexpected success when failure was expected
        (Ok(parsed_id), "fail", _) => {
            panic!(
                "Unexpected success for ID: {}. Expected outcome: fail, but got success with Assembly: '{:?}', Sequence ID: '{}', Ranges: {:?}, IDVersion: {:?}",
                id, parsed_id.assembly_id, parsed_id.sequence_id, parsed_id.ranges, parsed_id.inferred_version
            );
        }
        // Unexpected failure when pass was expected
        (Err(err), "pass", Some(expected_version)) => {
            panic!(
                "Unexpected failure for ID: {}. Expected outcome: pass, version: {:?}. Error: {}",
                id, expected_version, err
            );
        }
        // Catch-all for any unexpected outcome
        _ => panic!(
            "Unexpected outcome for ID: {}. Expected outcome: {}, version: {:?}, but got result: {:?}",
            id, expected_outcome, expected_version, result
        ),
    }
}
    #[test]
    fn test_parse_id() {
        let test_cases = vec![
            ("chr1:100-200_+", "pass", "", "chr1", vec![(100, 200, '+')], Some(IDVersion::V2), "chr1:100-200_+"),
            ("hg38:chr1:100-200_+", "pass", "hg38", "chr1", vec![(100, 200, '+')], Some(IDVersion::V2), "hg38:chr1:100-200_+"),
            ("hg38:chr1:100-200_-", "pass", "hg38", "chr1", vec![(100, 200, '-')], Some(IDVersion::V2), "hg38:chr1:100-200_-"),
            ("hg38:chr1:100-200_+:10-50_-:1-5_+", "pass", "hg38", "chr1", vec![(100, 200, '+'), (10, 50, '-'), (1, 5, '+')], Some(IDVersion::V2), "hg38:chr1:145-149_-"),
            ("hg38:chr1:200-100_-", "fail", "", "", vec![], None, ""),
            ("hg_38:chr+1:10-40_+", "pass", "hg_38", "chr+1", vec![(10, 40, '+')], Some(IDVersion::V2), "hg_38:chr+1:10-40_+"),
            // Add remaining test cases here in the same format
        ];

        for (id, outcome, assembly, seq_id, ranges, version, normalization) in test_cases {
            run_parse_id_test(id, outcome, assembly, seq_id, ranges, version, normalization);
        }
    }


}




