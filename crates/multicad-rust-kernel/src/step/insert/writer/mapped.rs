//! A single declared STEP map of the coupled local body; no world-point rewrite.
use super::*;
const OUTPUT: usize = 16 * 1024;
#[allow(clippy::too_many_arguments)]
pub fn write_mapped_candidate(
    topology: &TopologyPlan<'_, '_, '_>,
    origin: [f64; 3],
    u: [f64; 3],
    v: [f64; 3],
    maximum_work: usize,
    maximum_bytes: usize,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<String, Error> {
    let mut writer = Writer::new();
    writer.append_owned_insert_mapped_candidate(
        topology,
        origin,
        u,
        v,
        maximum_work,
        maximum_bytes,
        interrupt,
    )?;
    writer.finish().map_err(Error::from)
}
impl Writer {
    /// Append one local BRep with its ACTUAL representation map. Ratios encode
    /// the STEP normalized/projected basis, not the raw affine matrix [u,v,u×v].
    /// Native admission additionally requires the decoded member's material,
    /// mapping/source correspondence and independent mass proof.
    #[allow(clippy::too_many_arguments)]
    pub fn append_owned_insert_mapped_candidate(
        &mut self,
        topology: &TopologyPlan<'_, '_, '_>,
        origin: [f64; 3],
        u: [f64; 3],
        v: [f64; 3],
        maximum_work: usize,
        maximum_bytes: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<(), Error> {
        let n = pose(origin, u, v)?;
        self.append_typed(|file|{
            let (body,mut w,_memory)=encode_local(file,topology,maximum_work,maximum_bytes,interrupt,OUTPUT)?;
            // OUTPUT was reserved with the complete body before any write. All
            // following fixed records fit that bound; no dynamic geometry clone.
            w.charge(4096)?;
            let local_origin=point(file,[0.;3],&mut w)?;
            let source_axis=file.axis(local_origin,[0.,0.,1.],[1.,0.,0.]);
            let length=file.entity(format_args!("(LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.))"));
            let angle=file.entity(format_args!("(NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.))"));
            let solid_angle=file.entity(format_args!("(NAMED_UNIT(*) SI_UNIT($,.STERADIAN.) SOLID_ANGLE_UNIT())"));
            let context=file.entity(format_args!("(GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNIT_ASSIGNED_CONTEXT(({length},{angle},{solid_angle})) REPRESENTATION_CONTEXT('',''))"));
            let local=file.entity(format_args!("SHAPE_REPRESENTATION('',({body},{source_axis}),{context})"));
            let source=file.entity(format_args!("REPRESENTATION_MAP({source_axis},{local})"));
            let target_origin=point(file,origin,&mut w)?;
            let axis1=file.direction(u);let axis2=file.direction(v);let axis3=file.direction(n);
            let target=file.entity(format_args!("(CARTESIAN_TRANSFORMATION_OPERATOR({axis1},{axis2},{target_origin},1.)CARTESIAN_TRANSFORMATION_OPERATOR_3D({axis3})FUNCTIONALLY_DEFINED_TRANSFORMATION('','')GEOMETRIC_REPRESENTATION_ITEM()REPRESENTATION_ITEM(''))"));
            let mapped=file.entity(format_args!("MAPPED_ITEM('',{source},{target})"));w.charge(1)?;Ok(mapped)
        })
    }
}
fn enclosed(v: [f64; 3]) -> [I; 3] {
    v.map(|x| I::new(x.next_down(), x.next_up()))
}
fn dot(a: [I; 3], b: [I; 3]) -> I {
    a.into_iter()
        .zip(b)
        .fold(I::point(0.), |sum, (a, b)| sum + a * b)
}
pub(super) fn pose(origin: [f64; 3], u: [f64; 3], v: [f64; 3]) -> Result<[f64; 3], Error> {
    if origin.into_iter().chain(u).chain(v).any(|x| !x.is_finite()) {
        return Err(Error::Invalid);
    }
    // Enclosures cover the actual canonical 17-digit decimal ratios, not just
    // their binary64 caches. This domain bound guarantees a proper nonsingular
    // normalized STEP basis; it does NOT replace measured raw-frame transport.
    let [a, b] = [enclosed(u), enclosed(v)];
    let uu = dot(a, a) - I::point(1.);
    let vv = dot(b, b) - I::point(1.);
    let uv = dot(a, b);
    if [uu, vv, uv]
        .into_iter()
        .any(|x| !x.valid() || x.magnitude() > 1e-4)
    {
        return Err(Error::Unresolved);
    }
    let n = [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ];
    // The admitted Gram domain gives |u×v|² >= (1−δ)²−δ²
    // = 1−2δ > .99979 for the literal f64 domain threshold δ=1e−4.
    // Every input component is <1.001. Two products/subtraction,
    // input serialization and cross serialization together perturb each cross
    // component by <32eps (full ULP reserves also cover subnormals). Thus the
    // ACTUAL serialized determinant is >.99979−3*1.001*32eps >.9997.
    // A second determinant refusal is unreachable on this bounded domain.
    Ok(n)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::step_read::owned_insert::{CorrespondenceOptions, File, FluxOptions, Options};
    #[test]
    fn actual_mapped_producer_retains_one_closed_local_body_at_real_proper_poses() {
        for left in [false, true] {
            let body = crate::step_read::owned_insert::tests::graph(left, 4.);
            let d = DomainPlan::prepare(&body, DomainOptions::default(), None).unwrap();
            let s = SurfacePlan::prepare(&d, CurveOptions::default(), None).unwrap();
            let t = TopologyPlan::prepare(&s, 16_000_000, 64 * 1024 * 1024, None).unwrap();
            for (origin, u, v) in [
                ([0.; 3], [1., 0., 0.], [0., 1., 0.]),
                ([10., 20., 30.], [1., 0., 0.], [0., -1., 0.]),
                ([1e6, -2e6, 3e6], [0.8, 0.6, 0.], [0.36, -0.48, -0.8]),
            ] {
                let text =
                    write_mapped_candidate(&t, origin, u, v, 16_000_000, 64 * 1024 * 1024, None)
                        .unwrap();
                let file = File::parse(&text, Options::default(), None).unwrap();
                assert_eq!(file.represented_members().len(), 1);
                let member = file.represented_members()[0];
                assert!(member.mapping_link().is_some());
                let p = file
                    .member(member.item_entity(), Options::default(), None)
                    .unwrap();
                let b = p.support_bindings(Options::default(), None).unwrap();
                let material = b.material_bounds(Options::default(), None).unwrap();
                let mass = material
                    .occupied_volume(FluxOptions::default(), None)
                    .unwrap();
                assert!(mass.volume_mm3()[0] > 0.);
                let c = b
                    .source_correspondence(&t, CorrespondenceOptions::default(), None)
                    .unwrap();
                assert!(c.position_error_mm() < 1e-7);
                let map = p
                    .physical_mapping(
                        member.mapping_link().unwrap().operator_entity(),
                        Options::default(),
                        None,
                    )
                    .unwrap();
                assert_eq!(map.handedness().unwrap(), 1);
                let bounds = map.map_box_mm(material.bounds_mm()).unwrap();
                assert!(
                    bounds
                        .iter()
                        .all(|b| b[0].is_finite() && b[1].is_finite() && b[0] < b[1])
                );
            }
        }
    }
    #[test]
    fn mapped_pose_domain_and_atomic_compound_rollback_are_explicit() {
        let body = crate::step_read::owned_insert::tests::graph(false, 4.);
        let d = DomainPlan::prepare(&body, DomainOptions::default(), None).unwrap();
        let s = SurfacePlan::prepare(&d, CurveOptions::default(), None).unwrap();
        let t = TopologyPlan::prepare(&s, 16_000_000, 64 * 1024 * 1024, None).unwrap();
        assert!(matches!(
            write_mapped_candidate(
                &t,
                [f64::NAN, 0., 0.],
                [1., 0., 0.],
                [0., 1., 0.],
                16_000_000,
                64 * 1024 * 1024,
                None
            ),
            Err(Error::Invalid)
        ));
        let seed = || {
            let mut w = Writer::new();
            w.append_sphere(
                crate::sphere::Sphere::new(1.).unwrap(),
                [20., 0., 0.],
                [1., 0., 0.],
                [0., 1., 0.],
            )
            .unwrap();
            w
        };
        let expected = seed().finish().unwrap();
        for (origin, u, v, error) in [
            (
                [f64::NAN, 0., 0.],
                [1., 0., 0.],
                [0., 1., 0.],
                Error::Invalid,
            ),
            ([0.; 3], [1., 0., 0.], [1., 0., 0.], Error::Unresolved),
        ] {
            let mut w = seed();
            assert_eq!(
                w.append_owned_insert_mapped_candidate(
                    &t,
                    origin,
                    u,
                    v,
                    16_000_000,
                    64 * 1024 * 1024,
                    None
                ),
                Err(error)
            );
            assert_eq!(w.finish().unwrap(), expected);
        }
        let mut w = seed();
        assert_eq!(
            w.append_owned_insert_mapped_candidate(
                &t,
                [0.; 3],
                [1., 0., 0.],
                [0., 1., 0.],
                t.work_used() + 2000,
                64 * 1024 * 1024,
                None
            ),
            Err(Error::Limit)
        );
        assert_eq!(w.finish().unwrap(), expected);
        let mut w = seed();
        assert_eq!(
            w.append_owned_insert_mapped_candidate(
                &t,
                [0.; 3],
                [1., 0., 0.],
                [0., 1., 0.],
                16_000_000,
                64 * 1024 * 1024,
                Some(&|| true)
            ),
            Err(Error::Cancelled)
        );
        assert_eq!(w.finish().unwrap(), expected);
        // Both late failures occur after the coupled body was fully emitted.
        let (_, local_work, _) = encode_local(
            &mut StepFile::new(),
            &t,
            16_000_000,
            64 * 1024 * 1024,
            None,
            OUTPUT,
        )
        .unwrap();
        let mut late = seed();
        assert_eq!(
            late.append_owned_insert_mapped_candidate(
                &t,
                [0.; 3],
                [1., 0., 0.],
                [0., 1., 0.],
                local_work.used + 4096 + 256,
                64 * 1024 * 1024,
                None
            ),
            Err(Error::Limit)
        );
        assert_eq!(late.finish().unwrap(), expected);
        let polls = std::cell::Cell::new(0);
        let count = || {
            polls.set(polls.get() + 1);
            false
        };
        let mut complete = seed();
        complete
            .append_owned_insert_mapped_candidate(
                &t,
                [0.; 3],
                [1., 0., 0.],
                [0., 1., 0.],
                16_000_000,
                64 * 1024 * 1024,
                Some(&count),
            )
            .unwrap();
        let total = polls.get();
        let completed = complete.finish().unwrap();
        let polls = std::cell::Cell::new(0);
        let stop = || {
            polls.set(polls.get() + 1);
            polls.get() == total
        };
        let mut late = seed();
        assert_eq!(
            late.append_owned_insert_mapped_candidate(
                &t,
                [0.; 3],
                [1., 0., 0.],
                [0., 1., 0.],
                16_000_000,
                64 * 1024 * 1024,
                Some(&stop)
            ),
            Err(Error::Cancelled)
        );
        late.append_owned_insert_mapped_candidate(
            &t,
            [0.; 3],
            [1., 0., 0.],
            [0., 1., 0.],
            16_000_000,
            64 * 1024 * 1024,
            None,
        )
        .unwrap();
        assert_eq!(late.finish().unwrap(), completed);
    }
}
