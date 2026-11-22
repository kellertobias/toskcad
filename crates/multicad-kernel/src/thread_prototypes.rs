//! Source identity for reproducible geometry review reports.
use multicad_core::project_file::ContentHash;

/// Fingerprint of the Rust facade and its thread/source naming contracts.
/// Derived geometry assets and receipts are excluded.
pub fn source_fingerprint() -> String {
    let mut bytes = b"multicad.rust-thread-source/1".to_vec();
    for source in [
        crate::source_identity(),
        include_str!("threads.rs"),
        include_str!("geometry.rs"),
        include_str!("provenance.rs"),
        include_str!("lib.rs"),
        include_str!("provider_data.rs"),
    ] {
        bytes.extend_from_slice(&(source.len() as u64).to_le_bytes());
        bytes.extend_from_slice(source.as_bytes());
    }
    ContentHash::of(&bytes).to_string()
}
