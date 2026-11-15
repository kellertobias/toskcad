//! Fresh full-coefficient transfer into the very same coupled head/cavity chart.
//! The independent N and symbolic (z-topZ) are retained. No origin translation,
//! normal reconstruction, rigid projection or finite head/thread sewing follows.
use super::*;
pub struct CommonChartPeriodicThreadRepair<'a, 'input> {
    field: AnalyticThreadRepair<'a, 'input>,
}
impl<'a, 'input> CommonChartPeriodicThreadRepair<'a, 'input> {
    pub fn source(&self) -> &PeriodicThreadSource {
        self.field
            .periodic
            .as_ref()
            .expect("common periodic producer source")
    }
    pub fn chart(&self) -> CommonMouthChart {
        self.field.chart.expect("common periodic producer chart")
    }
    pub fn body(&self) -> &'a Body<'input> {
        self.field.body()
    }
    pub fn transfers(&self) -> &[AnalyticThreadTransfer] {
        self.field.transfers()
    }
    pub fn new_work_used(&self) -> usize {
        self.field.new_work_used()
    }
    pub fn inherited_head_work_used(&self) -> usize {
        self.field.inherited_head_work_used()
    }
    pub fn total_work_used(&self) -> usize {
        self.field.total_work_used()
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.field.retained_budget_bytes()
    }
}
impl<'h, 'q, 'r, 'a, 'input> ClosedRepairedHead<'h, 'q, 'r, 'a, 'input> {
    pub fn common_chart_periodic_thread_repair(
        &self,
        options: AnalyticThreadOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> std::result::Result<CommonChartPeriodicThreadRepair<'a, 'input>, PhillipsReadError> {
        self.build_thread_repair(options, interrupt, true, Some(self.chart()))
            .map(|field| CommonChartPeriodicThreadRepair { field })
    }
}
#[cfg(test)]
mod tests;
