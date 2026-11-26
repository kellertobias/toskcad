//! Versioned compact storage of exact Rust STEP source text: deflate of the
//! original bytes, in standard padded base64, behind a version prefix.
//! Expanding returns the original text byte for byte; it is lossless, not a
//! regenerated or canonicalised file.
use super::{DocumentError, invalid};

/// Version tag of the stored compact form.
pub const COMPACT_STEP_PREFIX: &str = "multicad.step.deflate.1:";
/// Largest STEP file the Rust readers accept, and so the largest text a
/// compact source may expand to.
pub const MAX_STEP_READ_BYTES: usize = 64 * 1024 * 1024;

const ALPHABET: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

fn encode(bytes: &[u8]) -> String {
    let mut out = String::with_capacity(bytes.len().div_ceil(3) * 4);
    for chunk in bytes.chunks(3) {
        let b = [
            chunk[0],
            chunk.get(1).copied().unwrap_or(0),
            chunk.get(2).copied().unwrap_or(0),
        ];
        let n = (u32::from(b[0]) << 16) | (u32::from(b[1]) << 8) | u32::from(b[2]);
        for i in 0..4 {
            if i <= chunk.len() {
                out.push(ALPHABET[(n >> (18 - 6 * i) & 63) as usize] as char);
            } else {
                out.push('=');
            }
        }
    }
    out
}
fn decode(text: &str) -> Option<Vec<u8>> {
    let bytes = text.as_bytes();
    if !bytes.len().is_multiple_of(4) {
        return None;
    }
    let value = |c: u8| ALPHABET.iter().position(|&a| a == c).map(|v| v as u32);
    let mut out = Vec::with_capacity(bytes.len() / 4 * 3);
    for (index, chunk) in bytes.chunks(4).enumerate() {
        let last = index + 1 == bytes.len() / 4;
        let padding = chunk.iter().rev().take_while(|&&c| c == b'=').count();
        if padding > 2 || (padding > 0 && !last) {
            return None;
        }
        let mut n = 0u32;
        for &c in &chunk[..4 - padding] {
            n = (n << 6) | value(c)?;
        }
        n <<= 6 * padding as u32;
        let decoded = [(n >> 16) as u8, (n >> 8) as u8, n as u8];
        // Canonical padding: the dropped low bits must be zero.
        if padding > 0 && decoded[3 - padding..].iter().any(|&b| b != 0) {
            return None;
        }
        out.extend_from_slice(&decoded[..3 - padding]);
    }
    Some(out)
}
/// Compact the exact original STEP text.
pub fn compact_step_text(text: &str) -> String {
    let deflated = miniz_oxide::deflate::compress_to_vec(text.as_bytes(), 9);
    let mut out = String::with_capacity(COMPACT_STEP_PREFIX.len() + deflated.len() / 3 * 4 + 4);
    out.push_str(COMPACT_STEP_PREFIX);
    out.push_str(&encode(&deflated));
    out
}
/// Whether `stored` has the versioned compact form (without expanding it).
pub fn is_compact_step_text(stored: &str) -> bool {
    stored.len() > COMPACT_STEP_PREFIX.len() && stored.starts_with(COMPACT_STEP_PREFIX)
}
/// Expand a stored compact source to its exact original text, refusing a
/// malformed form or one that would expand beyond the readers' own cap.
pub fn expand_step_text(stored: &str) -> Result<String, DocumentError> {
    let body = stored
        .strip_prefix(COMPACT_STEP_PREFIX)
        .ok_or_else(|| invalid("Compact STEP source has no supported version prefix"))?;
    let deflated =
        decode(body).ok_or_else(|| invalid("Compact STEP source is not valid base64"))?;
    let bytes = miniz_oxide::inflate::decompress_to_vec_with_limit(&deflated, MAX_STEP_READ_BYTES)
        .map_err(|_| {
            invalid(format!(
                "Compact STEP source is malformed or expands beyond {} MiB",
                MAX_STEP_READ_BYTES / (1024 * 1024)
            ))
        })?;
    String::from_utf8(bytes).map_err(|_| invalid("Compact STEP source is not UTF-8 text"))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn compact_step_text_round_trips_exact_bytes_and_refuses_malformed_forms() {
        for text in [
            "",
            "a",
            "ab",
            "abc",
            "ISO-10303-21;\n#1=CARTESIAN_POINT('',(1.E0,2.E0,3.E0));\n",
        ] {
            let stored = compact_step_text(text);
            assert!(stored.starts_with(COMPACT_STEP_PREFIX));
            assert_eq!(expand_step_text(&stored).unwrap(), text);
        }
        for bytes in [&b""[..], b"f", b"fo", b"foo", b"foob", b"fooba", b"foobar"] {
            assert_eq!(decode(&encode(bytes)).unwrap(), bytes);
        }
        assert_eq!(encode(b"foobar"), "Zm9vYmFy");
        assert_eq!(encode(b"fo"), "Zm8=");
        let good = compact_step_text("exact text");
        for bad in [
            "exact text".to_owned(),
            "multicad.step.deflate.2:AAAA".to_owned(),
            format!("{COMPACT_STEP_PREFIX}A"),
            format!("{COMPACT_STEP_PREFIX}Zm9=vYmF"),
            format!("{COMPACT_STEP_PREFIX}Zm9"),
            format!("{COMPACT_STEP_PREFIX}Zm9vYmFy"),
            good.replace('A', "$"),
        ] {
            assert!(expand_step_text(&bad).is_err(), "{bad}");
        }
        // Non-canonical padding bits are refused.
        assert!(decode("Zm9=").is_none());
        // A bomb expanding beyond the 64 MiB read cap is refused.
        let huge = vec![b' '; MAX_STEP_READ_BYTES + 1];
        let deflated = miniz_oxide::deflate::compress_to_vec(&huge, 9);
        let bomb = format!("{COMPACT_STEP_PREFIX}{}", encode(&deflated));
        assert!(expand_step_text(&bomb).is_err());
        let edge = String::from_utf8(vec![b' '; MAX_STEP_READ_BYTES]).unwrap();
        assert_eq!(
            expand_step_text(&compact_step_text(&edge)).unwrap().len(),
            MAX_STEP_READ_BYTES
        );
    }
}
