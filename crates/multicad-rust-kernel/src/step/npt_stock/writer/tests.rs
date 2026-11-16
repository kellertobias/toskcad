use super::*;
fn stopped_prefix(plan: &SurfacePlan<'_>, maximum_work: usize) -> (Result<(), Error>, usize) {
    let mut file = StepFile::new();
    let mut encoder = Encoder::new(&mut file, plan, maximum_work, 64 * 1024 * 1024, None).unwrap();
    let status = emit::body(&mut encoder, plan).map(|_| ());
    (status, encoder.file.next)
}
#[test]
fn real_record_transition_quotas_refuse_before_publishing_spline_support_edge_and_cap_records() {
    let source = super::super::tests::source(false);
    let options = CurveOptions::default();
    let plan = SurfacePlan::new(
        CurvePlan::new(&source, options, None).unwrap(),
        options.into(),
        None,
    )
    .unwrap();
    let file =
        write_source_candidate(&plan, options.maximum_work, options.maximum_bytes, None).unwrap();
    let id = |fragment: &str, occurrence: usize| -> usize {
        let record = file
            .as_str()
            .lines()
            .filter(|line| line.starts_with('#') && line.contains(fragment))
            .nth(occurrence)
            .unwrap();
        record[1..record.find('=').unwrap()].parse().unwrap()
    };
    let targets = [
        id("=CARTESIAN_POINT(", 0),
        id("=B_SPLINE_CURVE_WITH_KNOTS('',1,", 0),
        id("=DEFINITIONAL_REPRESENTATION(", 0),
        id("=AXIS2_PLACEMENT_3D(", 0),
        id("=SURFACE_CURVE(", 0),
        id("=EDGE_CURVE(", 0),
        id("=EDGE_LOOP(", 8),
    ];
    for target in targets {
        let mut bad = plan.work_used();
        let mut good = options.maximum_work;
        while good - bad > 1 {
            let middle = bad + (good - bad) / 2;
            let (_, next) = stopped_prefix(&plan, middle);
            if next <= target {
                bad = middle;
            } else {
                good = middle;
            }
        }
        let (status, next) = stopped_prefix(&plan, bad);
        assert_eq!(status, Err(Error::Limit));
        assert_eq!(
            next, target,
            "actual next record is the uncommitted selected boundary"
        );
        assert!(stopped_prefix(&plan, good).1 > target);
    }
}
