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
    exact_end: Option<crate::Rational>,
}
impl Source {
    pub fn treated(source: Cavity) -> Self {
        Self {
            thread: source.thread,
            entry: source.entry,
            tail: source.tail,
            intent: Some(source.intent),
            exact_end: source.exact_receiving_extent(),
        }
    }
    pub fn bare(thread: AuthoredThread) -> Self {
        Self {
            thread,
            entry: None,
            tail: None,
            intent: None,
            exact_end: thread.exact_axial_extent(),
        }
    }
    pub fn full_length_mm(&self) -> f64 {
        self.exact_end
            .map_or(self.thread.full_length_mm(), |x| x.as_f64())
    }
    pub fn exact_end(&self, w: &mut Work<'_>) -> Result<Q, Error> {
        if let Some(end) = self.exact_end {
            let (n, d) = end.ratio();
            Q::integer(n, w)?.div(&Q::integer(d, w)?, w)
        } else {
            Q::float(self.full_length_mm(), w)
        }
    }
    pub fn chart(&self) -> ThreadChart<'_> {
        ThreadChart {
            thread: &self.thread,
            entry: self.entry,
            top: self.tail.map_or(self.full_length_mm(), |t| t.station),
            exact_top: if self.tail.is_none() {
                self.exact_end
            } else {
                None
            },
        }
    }
}
