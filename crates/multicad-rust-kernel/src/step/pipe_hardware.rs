//! Authored finite-head STEP proposals and independent qualification.
//! No Model/Native capability is granted by this module.
pub(crate) mod numbers;
mod writer;
pub use writer::{CandidateFile, Options, Receipt, propose};
mod qualify;
pub use qualify::{QualificationOptions, QualificationReceipt, QualifiedStep, qualify_actual_file};
#[cfg(test)]
mod tests;

#[cfg(test)]
pub(crate) mod test_fixture;
