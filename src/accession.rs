use std::sync::OnceLock;

use rayon::prelude::*;
use regex::Regex;

fn accession_regex() -> &'static Regex {
    static RE: OnceLock<Regex> = OnceLock::new();
    RE.get_or_init(|| {
        Regex::new(concat!(
            r"^(?:",
            r"[A-Z]{2}_[A-Z]{2}[0-9]{5,}|",
            r"[A-Z]{2}_[A-Z]{1,6}[0-9]{5,}(?:[A-Z]{0,2})?|",
            r"[A-Z]{1,4}_?[0-9]{5,}|",
            r"[A-Z]{4,6}[0-9]{8,}(?:[A-Z]{0,2})?|",
            r"[A-Z]{3}[0-9]{5}|",
            r"[A-Z][0-9][A-Z0-9]{8}|",
            r"[A-Z][0-9][A-Z0-9]{3}[0-9]",
            r")(?:(?:\.)([0-9]+))?"
        ))
        .expect("accession regex is valid")
    })
}

/// Extract the first accession from a string.
///
/// Returns `"NA"` when no accession is found, matching Python `taxutils`.
pub fn parse_accession(text: &str, version: bool) -> String {
    let upper = text.to_ascii_uppercase();
    for (start, _) in upper.char_indices() {
        if start > 0 {
            let previous = upper.as_bytes()[start - 1];
            if previous.is_ascii_alphanumeric() || previous == b'_' {
                continue;
            }
        }
        let Some(captures) = accession_regex().captures(&upper[start..]) else {
            continue;
        };
        let whole = captures.get(0).expect("whole match");
        let end = start + whole.end();
        if upper
            .as_bytes()
            .get(end)
            .is_some_and(|c| c.is_ascii_alphanumeric() || *c == b'_')
        {
            continue;
        }
        let mut accession = whole.as_str().to_owned();
        if !version && captures.get(1).is_some() {
            accession.truncate(accession.rfind('.').expect("version separator"));
        }
        return accession;
    }
    "NA".to_owned()
}

/// Batch form of [`parse_accession`].
pub fn parse_accessions<S: AsRef<str> + Sync>(strings: &[S], version: bool) -> Vec<String> {
    strings
        .par_iter()
        .map(|text| parse_accession(text.as_ref(), version))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_python_examples() {
        assert_eq!(
            parse_accession(">NC_045512.2 SARS-CoV-2", true),
            "NC_045512.2"
        );
        assert_eq!(
            parse_accession(">kraken:taxid|2886930|NC_001422.1 phiX", true),
            "NC_001422.1"
        );
        assert_eq!(parse_accession(">nc_045512.2", false), "NC_045512");
        assert_eq!(parse_accession("not_an_accession", true), "NA");
    }

    #[test]
    fn respects_identifier_boundaries() {
        assert_eq!(parse_accession("FOONC_045512.2", true), "NA");
        assert_eq!(parse_accession("NC_045512XYZ", true), "NA");
    }
}
