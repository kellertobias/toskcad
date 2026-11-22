//! Actual source oracles and bounded-work controls for production root caps.
use super::bound_root::{
    Candidate as Probe, EndpointArithmetic, decode, eval, parameter_envelope, source_gate,
    source_widths,
};
pub(super) use super::jets::field_bounds;
use super::*;
use crate::authored_hanger_bolt::{BoundaryOptions, Hanger};
fn run(curve: &CapCurvePiece<'_, '_>) -> Result<Probe, Error> {
    run_with_limit(curve, 16_000_000)
}
pub(super) fn run_with_limit(curve: &CapCurvePiece<'_, '_>, limit: usize) -> Result<Probe, Error> {
    let proof = run_bound(curve, None, limit, None, false)?;
    source_gate(&proof)?;
    Ok(proof)
}
pub(super) fn run_to_vertices(
    curve: &CapCurvePiece<'_, '_>,
    targets: [&[EncodedScalar; 3]; 2],
    limit: usize,
) -> Result<Probe, Error> {
    let proof = run_bound(curve, Some(targets), limit, None, false)?;
    source_gate(&proof)?;
    Ok(proof)
}
fn run_bound(
    curve: &CapCurvePiece<'_, '_>,
    targets: Option<[&[EncodedScalar; 3]; 2]>,
    limit: usize,
    interrupt: Option<&dyn Fn() -> bool>,
    exact_endpoints: bool,
) -> Result<Probe, Error> {
    let mut w = Work::new(limit, interrupt)?;
    w.charge(curve.field.work_used() + curve.work_used)?;
    let mut m = Memory::new(
        curve.field.retained_bytes() + curve.retained_bytes,
        64 * 1024 * 1024,
    )?;
    m.add(SCRATCH)?;
    super::bound_root::bind(
        curve,
        targets,
        &mut m,
        &mut w,
        if exact_endpoints {
            EndpointArithmetic::Exact
        } else {
            EndpointArithmetic::Adaptive
        },
    )
}
#[test]
fn actual_root_full_basis_and_same_tensor_pcurve_lift_pass_both_hands_nonbinary() {
    for (left, nb, period, end, piece) in [
        (false, false, -1, End::Start, 3),
        (true, true, 25, End::End, 0),
    ] {
        let s = Hanger::new(crate::authored_hanger_bolt::tests::fixture(left, nb), None).unwrap();
        let b = s.boundary(BoundaryOptions::default(), None).unwrap();
        let f = EncodedField::prepare(&b, CurveOptions::default(), None).unwrap();
        let i = f
            .cells()
            .iter()
            .position(|c| {
                c.source_key.total_period == period
                    && matches!(
                        c.source_key.face,
                        FaceKey::MetricThread {
                            band: ThreadBandRole::Root,
                            ..
                        }
                    )
            })
            .unwrap();
        let curve = f
            .cap_curve_piece(i, end, piece, DomainOptions::default(), None)
            .unwrap()
            .unwrap();
        let proof = run(&curve).unwrap();
        independent_oracle(&curve, &proof);
        let targets = std::array::from_fn(|e| {
            let mut memory = Memory::new(0, 64 * 1024 * 1024).unwrap();
            let mut work = Work::new(16_000_000, None).unwrap();
            let u = &proof.parameter[e];
            let xyz = eval(&decode(curve.controls(), &mut work).unwrap(), u, &mut work).unwrap();
            std::array::from_fn(|axis| {
                EncodedScalar::retain(
                    &xyz[axis].div(&xyz[3], &mut work).unwrap(),
                    &mut memory,
                    &mut work,
                )
                .unwrap()
            })
        });
        let bound = run_to_vertices(&curve, targets.each_ref(), 16_000_000).unwrap();
        assert!(bound.vertex_error < 1e-7);
        eprintln!(
            "actual full-basis+UV left={left} nb={nb} cumulative_work={} position={} derivative={}",
            proof.work, proof.position, proof.derivative
        );
    }
}

#[test]
fn old_meridian_u_chart_retains_actual_field_lift_and_real_source_c1_refusal() {
    for (left, nb) in [(false, false), (true, true)] {
        let source =
            Hanger::new(crate::authored_hanger_bolt::tests::fixture(left, nb), None).unwrap();
        let boundary = source.boundary(BoundaryOptions::default(), None).unwrap();
        let field = EncodedField::prepare(&boundary, CurveOptions::default(), None).unwrap();
        let mut count = 0;
        let mut subtotal = field.work_used();
        let mut worst = [0_f64; 2];
        let mut stages = [0usize; 7];
        let mut vertex = 0_f64;
        let mut source_refusals = 0;
        for index in 0..field.cells().len() {
            for end in [End::Start, End::End] {
                for piece in 0..field.angular.coefficients().len() {
                    let Some(curve) = field
                        .cap_curve_piece(index, end, piece, DomainOptions::default(), None)
                        .unwrap()
                    else {
                        continue;
                    };
                    // The historical chart is a paired diagnostic, not the
                    // new Source-qualified carrier route. Preserve its real
                    // tiny-span C1 refusal while checking the field proof.
                    let proof = run_bound(&curve, None, 16_000_000, None, false).unwrap();
                    if source_gate(&proof).is_err() {
                        source_refusals += 1;
                        assert!(nb);
                        assert!(matches!(
                            field.cells[index].source_key.face,
                            FaceKey::WoodThread {
                                band: ThreadBandRole::Crest,
                                ..
                            }
                        ));
                        eprintln!(
                            "old-u source C1 refusal left={left} nb={nb} cell={index} end={end:?} piece={piece} originalD1={}",
                            proof.original_source_derivative
                        );
                    }
                    subtotal += proof.work - field.work_used();
                    for (sum, value) in stages.iter_mut().zip(proof.stages) {
                        *sum += value;
                    }
                    vertex = vertex.max(proof.vertex_error);
                    worst[0] = worst[0].max(proof.position);
                    worst[1] = worst[1].max(proof.derivative);
                    count += 1;
                }
            }
        }
        eprintln!(
            "old-u field-only source_refusals={source_refusals} caps-only positive-query subtotal left={left} nb={nb} count={count} work={subtotal} worst={worst:?} vertex={vertex} stages={stages:?}; excludes selection/field-surfaces/closedfile"
        );
        assert!(count >= 44);
        assert_eq!(source_refusals > 0, nb);
    }
}

fn independent_oracle(curve: &CapCurvePiece<'_, '_>, proof: &Probe) {
    use super::restrict::tests::jet;
    use crate::step::hanger::encoded::cap::tests::{bezier, number};
    use num_rational::BigRational as Q;
    use num_traits::Signed;
    let homogeneous = |data: &super::literal::CartesianCap| {
        data.poles
            .iter()
            .zip(&data.weights)
            .map(|(p, w)| {
                let p = p.each_ref().map(number);
                let w = number(w);
                [&p[0] * &w, &p[1] * &w, &p[2] * &w, w]
            })
            .collect::<Vec<_>>()
    };
    let xyz = homogeneous(&proof.xyz);
    let uv = homogeneous(&proof.uv);
    let f = curve.field;
    let row = &f.meridians[f.cells[curve.source_cell].meridian];
    let controls = (0..=row.degree())
        .map(|i| {
            (0..=10)
                .map(|j| {
                    f.control(
                        curve.source_cell,
                        i,
                        curve.angular_piece,
                        j,
                        DomainOptions::default(),
                        None,
                    )
                    .unwrap()
                    .homogeneous
                    .each_ref()
                    .map(number)
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let ka = Q::from_float(f.angular.knots()[curve.angular_piece]).unwrap();
    let kw = Q::from_float(f.angular.knots()[curve.angular_piece + 1]).unwrap() - &ka;
    let span = &proof.parameter[1].0 - &proof.parameter[0].0;
    for k in 0..=4 {
        let t = Q::new(k.into(), 4.into());
        let u3 = &proof.parameter[0].0 + &span * &t;
        let (point, derivative) = jet(xyz.clone(), &u3);
        let derivative = derivative.map(|d| d * &span);
        let (foot, foot_derivative) = jet(uv.clone(), &t);
        let v = (&foot[1] - &ka) / &kw;
        let hrows = controls
            .iter()
            .map(|r| bezier(r.clone(), &v))
            .collect::<Vec<_>>();
        let h = bezier(hrows.clone(), &foot[0]);
        let hu = bezier(
            hrows
                .windows(2)
                .map(|p| {
                    std::array::from_fn::<_, 4, _>(|axis| {
                        (&p[1][axis] - &p[0][axis]) * Q::from_integer(row.degree().into())
                    })
                })
                .collect(),
            &foot[0],
        );
        let hvrows = controls
            .iter()
            .map(|r| {
                bezier(
                    r.windows(2)
                        .map(|p| {
                            std::array::from_fn::<_, 4, _>(|axis| {
                                (&p[1][axis] - &p[0][axis]) * Q::from_integer(10.into()) / &kw
                            })
                        })
                        .collect(),
                    &v,
                )
            })
            .collect::<Vec<_>>();
        let hv = bezier(hvrows, &foot[0]);
        let surface = std::array::from_fn::<_, 3, _>(|axis| &h[axis] / &h[3]);
        let surface_derivative = std::array::from_fn::<_, 3, _>(|axis| {
            ((&hu[axis] * &h[3] - &h[axis] * &hu[3]) * &foot_derivative[0]
                + (&hv[axis] * &h[3] - &h[axis] * &hv[3]) * &foot_derivative[1])
                / (&h[3] * &h[3])
        });
        let gap = (0..3).map(|i| (&point[i] - &surface[i]).abs()).sum::<Q>();
        let tangent = (0..3)
            .map(|i| (&derivative[i] - &surface_derivative[i]).abs())
            .sum::<Q>();
        assert!(gap <= Q::from_float(proof.position).unwrap());
        assert!(tangent <= Q::from_float(proof.derivative).unwrap());
    }
}

#[test]
fn same_actual_root_exact_endpoint_fallback_and_every_earned_checkpoint_remain_bounded() {
    use std::cell::Cell;
    for x in [
        [f64::NAN, 0., 0.],
        [-1., 0., 0.],
        [0., 2e-6, 0.],
        [0., 0., 2e-6],
    ] {
        assert!(matches!(source_widths(x), Err(Error::Unresolved)));
    }
    source_widths([0., 0., 0.]).unwrap();
    let s = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(false, false),
        None,
    )
    .unwrap();
    let b = s.boundary(BoundaryOptions::default(), None).unwrap();
    let f = EncodedField::prepare(&b, CurveOptions::default(), None).unwrap();
    let index = f
        .cells()
        .iter()
        .position(|c| {
            c.source_key.total_period == -1
                && matches!(
                    c.source_key.face,
                    FaceKey::MetricThread {
                        band: ThreadBandRole::Root,
                        ..
                    }
                )
        })
        .unwrap();
    let curve = f
        .cap_curve_piece(index, End::Start, 3, DomainOptions::default(), None)
        .unwrap()
        .unwrap();
    let interval = run(&curve).unwrap();
    let exact = run_bound(&curve, None, 16_000_000, None, true).unwrap();
    source_gate(&exact).unwrap();
    independent_oracle(&curve, &exact);
    let mut memory = Memory::new(0, 64 * 1024 * 1024).unwrap();
    let mut work = Work::new(16_000_000, None).unwrap();
    let h = decode(curve.controls(), &mut work).unwrap();
    let targets = interval.parameter.each_ref().map(|u| {
        let p = eval(&h, u, &mut work).unwrap();
        std::array::from_fn::<_, 3, _>(|a| {
            EncodedScalar::retain(&p[a].div(&p[3], &mut work).unwrap(), &mut memory, &mut work)
                .unwrap()
        })
    });
    let exact_target = run_bound(&curve, Some(targets.each_ref()), 16_000_000, None, true).unwrap();
    source_gate(&exact_target).unwrap();
    // A real source-root envelope must also cover an explicitly coarser
    // parameter proposal. This is a numeric hull check, not curve admission.
    let half_root = f
        .cap_curve_piece(index, End::Start, 4, DomainOptions::default(), None)
        .unwrap()
        .unwrap();
    let proposed = [N::integer(0), N::integer(1)];
    let envelope = parameter_envelope(half_root.parameter_bounds, &proposed, &mut work).unwrap();
    assert_eq!(envelope[0].0, proposed[0].0);
    assert_eq!(envelope[1].0, proposed[1].0);
    assert!(
        exact.work > interval.work,
        "attempted interval proposal and exact fallback are both charged"
    );
    let calls = Cell::new(0);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    run_bound(&curve, None, 16_000_000, Some(&count), false).unwrap();
    for stop in 1..=calls.get() {
        let n = Cell::new(0);
        let cancel = || {
            n.set(n.get() + 1);
            n.get() >= stop
        };
        assert_eq!(
            run_bound(&curve, None, 16_000_000, Some(&cancel), false).err(),
            Some(Error::Cancelled),
            "stop={stop}"
        );
    }
    eprintln!(
        "actual root carrier checkpoint callbacks={} intervalwork={} exactwork={}",
        calls.get(),
        interval.work,
        exact.work
    );
}
