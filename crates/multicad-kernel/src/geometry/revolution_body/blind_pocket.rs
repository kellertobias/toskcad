//! Retained source-certified coaxial blind pockets, never inferred from a lathe.
use super::*;
pub use multicad_rust_kernel::cylindrical_blind_pocket::Spec as CylindricalBlindPocketSpec;
use multicad_rust_kernel::cylindrical_blind_pocket::{Frame, Placement, Pocket, Role};
impl Solid {
    pub fn rust_cylindrical_stock_available(&self) -> bool {
        self.axial_cylinder_frame().is_some()
    }

    pub fn rust_cut_cylindrical_blind_pocket(
        &self,
        spec: CylindricalBlindPocketSpec,
        entry: PlaneFrame,
        feature: EntityId,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, GeometryError> {
        let check = || {
            if interrupt.is_some_and(|f| f()) {
                Err(error(
                    "cylindrical blind pocket",
                    GeometryErrorKind::Cancelled,
                    "Receiving cut superseded",
                ))
            } else {
                Ok(())
            }
        };
        check()?;
        policy.validate()?;
        entry.checked(policy)?;
        if self.policy != policy || policy.boolean_fuzzy_mm != 0. || feature.is_nil() {
            return Err(input(
                "Cylindrical blind pockets require a distinct identity and the stock's exact policy",
            ));
        }
        let Some((stock, stock_frame)) = self.axial_cylinder_frame() else {
            return Err(error(
                OPERATION,
                GeometryErrorKind::Unsupported,
                "Blind cylinder pockets require one unmachined actual stock cylinder",
            ));
        };
        if stock.feature() == feature.uuid_bits() {
            return Err(input("Receiving and stock source identities must differ"));
        }
        let raw = |f: PlaneFrame| Frame {
            origin_mm: f.origin,
            axes: [f.u, f.v, f.normal()],
        };
        let placement = if self.axial_cylinder().is_some() {
            Placement::Axial
        } else {
            Placement::Framed(raw(stock_frame))
        };
        let source = Pocket::new(
            *stock,
            placement,
            raw(entry),
            spec,
            policy.linear_tolerance_mm,
            policy.minimum_feature_mm,
            interrupt,
        )
        .map_err(source_failure)?;
        let raw = source.frame();
        let frame = PlaneFrame {
            origin: raw.origin_mm,
            u: raw.axes[0],
            v: raw.axes[1],
        };
        // The stock's constructors and transformed() already proved this
        // exact radius/radius/height metric domain. Framed retains its literal
        // columns; the axial-WCS case uses exact cardinal material columns.
        let geometry = source.geometry().clone();
        let original = self.topology_snapshot(feature, 0)?;
        let cap = |end| {
            original.elements.iter().find(|e| {
                e.key.kind == multicad_core::references::TopologyKind::Face
                    && e.key.ordinal == u32::from(end)
            })
        };
        let outer = original.elements.iter().find(|e| {
            e.key.kind == multicad_core::references::TopologyKind::Face && e.key.ordinal == 2
        });
        let mut segments = vec![];
        for (&role, &authored_start) in source.roles().iter().zip(&geometry.rings()[0]) {
            let (name, ancestor) = match role {
                Role::Back => ("stock.back", cap(!source.end())),
                Role::Outer => ("stock.outer", outer),
                Role::Mouth => ("entry", cap(source.end())),
                Role::SeatWall => ("seat.wall", None),
                Role::Shoulder => ("seat.shoulder", None),
                Role::Bore => ("bore", None),
                Role::Floor => ("floor", None),
                Role::Axis => ("axis", None),
            };
            let inherited_frames = ancestor
                .into_iter()
                .flat_map(|e| &e.frames)
                .map(|a| {
                    Ok(AnchoredFrame {
                        origin: a.origin.clone(),
                        frame: local_attachment(frame, a.frame, policy)?,
                    })
                })
                .collect::<Result<Vec<_>, GeometryError>>()?;
            segments.push(MeridianSegment {
                origin: Origin {
                    source: feature,
                    role: format!("connection.receiver/{name}"),
                },
                direct_identity: true,
                inherited: ancestor
                    .into_iter()
                    .flat_map(|e| e.origins.iter().cloned())
                    .collect(),
                inherited_frames,
                authored_start,
            });
        }
        let mut retained = vec![];
        let mut ordinal = 0;
        for &[r, z] in &geometry.rings()[0] {
            if r > 0. {
                if r == stock.circle().radius().as_f64() {
                    let end = z > 0.;
                    if let Some(ancestor) = cap(end)
                        && let Some(edge) = original.elements.iter().find(|e| {
                            e.key.kind == multicad_core::references::TopologyKind::Edge
                                && ancestor.origins.iter().any(|o| e.boundary.contains(o))
                        })
                    {
                        retained.push((ordinal, edge.boundary.clone()));
                    }
                }
                ordinal += 1;
            }
        }
        let body = RevolutionBody {
            rings: topology::ring_edges(&geometry),
            geometry,
            feature,
            sketch: PlaneFrame {
                v: [0., 0., 1.],
                ..Default::default()
            },
            segments: vec![segments],
            reverse_rims: false,
            retained_rim_boundaries: retained,
            cylindrical_pocket: Some(Rc::new(source)),
        };
        let bounds = body.bounds(frame);
        policy.point(bounds.min)?;
        policy.point(bounds.max)?;
        check()?;
        Ok(Self {
            body: SolidBody::Revolved(Rc::new(body), frame),
            policy,
            provenance: vec![],
            rust_origins: vec![],
            retained_placement: None,
        })
    }
    pub fn rust_cylindrical_blind_pocket_geometry(&self) -> Option<(&Pocket, PlaneFrame)> {
        let SolidBody::Revolved(body, frame) = &self.body else {
            return None;
        };
        Some((body.cylindrical_pocket.as_deref()?, *frame))
    }
}
impl RevolutionBody {
    pub(crate) fn has_cylindrical_blind_pocket(&self) -> bool {
        self.cylindrical_pocket.is_some()
    }
}
#[cfg(test)]
mod tests;

fn source_failure(cause: multicad_rust_kernel::tapered_thread::Error) -> GeometryError {
    use multicad_rust_kernel::tapered_thread::Error;
    let kind = match cause {
        Error::Invalid => GeometryErrorKind::Input,
        Error::Cancelled => GeometryErrorKind::Cancelled,
        Error::Limit => GeometryErrorKind::Limit,
        _ => GeometryErrorKind::Unsupported,
    };
    error(
        "cylindrical blind pocket",
        kind,
        format!("Actual cylinder pocket certificate refused: {cause:?}"),
    )
}
