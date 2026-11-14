//! Factual carrier chart shared by genuine treated and bare source owners.
//! This is not a no-op Cavity; its optional data is copied only from a real
//! treatment, and the public Boundary separately retains the original type.
use super::*;
#[derive(Clone, Copy, PartialEq)]
pub(super) struct Source {
    pub thread: AuthoredThread,
    pub entry: Option<Entry>,
    pub tail: Option<Tail>,
    // Preserve the complete actual treatment identity for stock-material checks.
    intent: Option<Intent>,
}
impl Source {
    pub fn treated(source: Cavity) -> Self {
        Self {
            thread: source.thread,
            entry: source.entry,
            tail: source.tail,
            intent: Some(source.intent),
        }
    }
    pub fn bare(thread: AuthoredThread) -> Self {
        Self {
            thread,
            entry: None,
            tail: None,
            intent: None,
        }
    }
    pub fn full_length_mm(&self) -> f64 {
        self.thread.full_length_mm()
    }
    pub fn chart(&self) -> ThreadChart<'_> {
        ThreadChart {
            thread: &self.thread,
            entry: self.entry,
            top: self.tail.map_or(self.full_length_mm(), |t| t.station),
        }
    }
}
