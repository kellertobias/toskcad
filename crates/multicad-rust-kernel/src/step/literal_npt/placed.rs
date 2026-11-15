//! Flatten ONLY the output operator; the retained original geometry/map stay immutable.
use super::*;
struct Local<'a> {
    source: &'a OwnedMaterial,
    indices: Vec<usize>,
}
impl RetainedRecords for Local<'_> {
    fn record_count(&self) -> usize {
        self.indices.len()
    }
    fn record_index(&self, id: u32) -> Option<usize> {
        self.indices
            .binary_search(&self.source.record_index(id)?)
            .ok()
    }
    fn component(&self, record: usize, index: usize) -> Option<&str> {
        self.source.component(*self.indices.get(record)?, index)
    }
}
fn columns(
    source: &OwnedMaterial,
    pose: &read::LiteralOccurrencePlacement,
) -> ([[I; 3]; 3], [I; 3]) {
    let old = source
        .linear_columns_bounds()
        .map(|c| c.map(|b| I::new(b[0], b[1])));
    let matrix = pose
        .linear_rows_bounds()
        .map(|r| r.map(|b| I::new(b[0], b[1])));
    let origin = source.origin_mm().map(|b| I::new(b[0], b[1]));
    let shift = pose.translation_bounds_mm().map(|b| I::new(b[0], b[1]));
    (
        std::array::from_fn(|j| std::array::from_fn(|i| dot(matrix[i], old[j]))),
        std::array::from_fn(|i| dot(matrix[i], origin) + shift[i]),
    )
}
fn dot(a: [I; 3], b: [I; 3]) -> I {
    a.into_iter()
        .zip(b)
        .fold(I::point(0.), |s, (a, b)| s + a * b)
}
fn cross(a: [I; 3], b: [I; 3]) -> [I; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn l1(a: [I; 3]) -> I {
    a.into_iter()
        .fold(I::point(0.), |s, v| s + I::point(v.magnitude()))
}
fn norm(a: [I; 3]) -> I {
    a.into_iter()
        .fold(I::point(0.), |s, v| s + v.square())
        .sqrt()
}
fn entity(file: &mut StepFile, trace: &mut Trace, body: fmt::Arguments<'_>) -> Entity {
    let prefix = file.data.len();
    let id = file.entity(body);
    trace.nodes.push(Node {
        id: u32::try_from(id.0).expect("bounded STEP entity"),
        body: prefix + digits(id.0) + 2..file.data.len() - 2,
    });
    id
}
pub(super) fn encode(
    source: &OwnedMaterial,
    pose: &read::LiteralOccurrencePlacement,
    file: &mut StepFile,
    w: &mut Ledger<'_>,
) -> Result<(Trace, usize)> {
    let before = w.used;
    let count = source.record_count();
    // Initial selection initialization, complete index scan and bounded
    // selected-index publication are actual new preparation work.
    w.step(count.checked_mul(3).ok_or(Error::Limit)?)?;
    let requested = count
        .checked_mul(1 + 2 * size_of::<usize>())
        .ok_or(Error::Limit)?;
    w.observe(
        source.retained_budget_bytes()
            + pose.retained_budget_bytes()
            + file.data.capacity()
            + requested
            + FIXED,
    )?;
    let mut selected = Vec::new();
    selected
        .try_reserve_exact(count)
        .map_err(|_| Error::Limit)?;
    selected.resize(count, false);
    let root = source
        .record_index(source.member().local_solid_entity())
        .expect("retained local root");
    let mut pending = Vec::new();
    pending.try_reserve_exact(count).map_err(|_| Error::Limit)?;
    pending.push(root);
    selected[root] = true;
    let mut indices = Vec::new();
    indices
        .try_reserve_exact(source.record_count())
        .map_err(|_| Error::Limit)?;
    w.observe(
        source.retained_budget_bytes()
            + pose.retained_budget_bytes()
            + file.data.capacity()
            + selected.capacity()
            + pending.capacity() * size_of::<usize>()
            + indices.capacity() * size_of::<usize>()
            + FIXED,
    )?;
    while let Some(index) = pending.pop() {
        w.step(32)?;
        let mut component = 0;
        while let Some(mut value) = source.component(index, component) {
            w.step(value.len())?;
            while let Some(value) = token(&mut value) {
                if let Token::Reference(id) = value {
                    w.step(32)?;
                    let next = source
                        .record_index(id)
                        .expect("retained closed dependency graph");
                    if !selected[next] {
                        selected[next] = true;
                        pending.push(next);
                    }
                }
            }
            component += 1;
        }
    }
    indices.extend(
        selected
            .iter()
            .enumerate()
            .filter_map(|(i, v)| v.then_some(i)),
    );
    let local = Local { source, indices };
    let resident = source.retained_budget_bytes()
        + pose.retained_budget_bytes()
        + selected.capacity()
        + pending.capacity() * size_of::<usize>()
        + local.indices.capacity() * size_of::<usize>();
    w.observe(resident + file.data.capacity() + FIXED)?;
    let graph_work = w.used - before;
    let mut trace = super::encode(
        &local,
        source.member().local_solid_entity(),
        file,
        resident,
        w,
    )?;
    let body = Entity(trace.member as usize);
    w.observe(
        resident
            + file.data.capacity().max(file.data.len() + OUTPUT)
            + trace.nodes.capacity().max(trace.nodes.len() + 24) * size_of::<Node>()
            + FIXED,
    )?;
    trace
        .nodes
        .try_reserve_exact(24)
        .map_err(|_| Error::Limit)?;
    file.data
        .try_reserve_exact(OUTPUT)
        .map_err(|_| Error::Limit)?;
    w.step(8192)?;
    let (columns, origin) = columns(source, pose);
    let scale = norm(columns[0]).midpoint();
    if !scale.is_finite() || scale <= 0. {
        return Err(Error::Unsupported);
    }
    let axes = columns.map(|c| c.map(|v| v.midpoint() / scale));
    let zero = entity(
        file,
        &mut trace,
        format_args!("CARTESIAN_POINT('',(0.,0.,0.))"),
    );
    let z = entity(file, &mut trace, format_args!("DIRECTION('',(0.,0.,1.))"));
    let x = entity(file, &mut trace, format_args!("DIRECTION('',(1.,0.,0.))"));
    let axis = entity(
        file,
        &mut trace,
        format_args!("AXIS2_PLACEMENT_3D('',{zero},{z},{x})"),
    );
    let length = entity(
        file,
        &mut trace,
        format_args!("(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.))"),
    );
    let angle = entity(
        file,
        &mut trace,
        format_args!("(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.))"),
    );
    let solid_angle = entity(
        file,
        &mut trace,
        format_args!("(NAMED_UNIT(*)SI_UNIT($,.STERADIAN.)SOLID_ANGLE_UNIT())"),
    );
    let context = entity(
        file,
        &mut trace,
        format_args!(
            "(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNIT_ASSIGNED_CONTEXT(({length},{angle},{solid_angle}))REPRESENTATION_CONTEXT('',''))"
        ),
    );
    let rep = entity(
        file,
        &mut trace,
        format_args!("SHAPE_REPRESENTATION('',({body},{axis}),{context})"),
    );
    let map = entity(
        file,
        &mut trace,
        format_args!("REPRESENTATION_MAP({axis},{rep})"),
    );
    let target = entity(
        file,
        &mut trace,
        format_args!(
            "CARTESIAN_POINT('',{})",
            Coordinates(origin.map(I::midpoint))
        ),
    );
    let a = entity(
        file,
        &mut trace,
        format_args!("DIRECTION('',{})", Coordinates(axes[0])),
    );
    let b = entity(
        file,
        &mut trace,
        format_args!("DIRECTION('',{})", Coordinates(axes[1])),
    );
    let c = entity(
        file,
        &mut trace,
        format_args!("DIRECTION('',{})", Coordinates(axes[2])),
    );
    let operator = entity(
        file,
        &mut trace,
        format_args!(
            "(CARTESIAN_TRANSFORMATION_OPERATOR({a},{b},{target},{scale:.16E})CARTESIAN_TRANSFORMATION_OPERATOR_3D({c})FUNCTIONALLY_DEFINED_TRANSFORMATION('','')GEOMETRIC_REPRESENTATION_ITEM()REPRESENTATION_ITEM(''))"
        ),
    );
    let item = entity(
        file,
        &mut trace,
        format_args!("MAPPED_ITEM('',{map},{operator})"),
    );
    trace.member = u32::try_from(item.0).expect("bounded STEP mapped item");
    w.observe(
        resident + file.data.capacity() + trace.nodes.capacity() * size_of::<Node>() + FIXED,
    )?;
    if w.used - before - graph_work > 16_000_000 {
        return Err(Error::Limit);
    }
    Ok((trace, graph_work))
}
pub(super) fn correspond(
    source: &OwnedMaterial,
    placement: &read::LiteralOccurrencePlacement,
    decoded: &read::Placement<'_, '_, '_>,
    options: Options,
    w: &mut Ledger<'_>,
) -> Result<(f64, f64, f64)> {
    w.step(8192)?;
    let (wanted, origin) = columns(source, placement);
    let actual = decoded
        .linear_columns_bounds()
        .map(|r| r.map(|b| I::new(b[0], b[1])));
    let actual_origin = decoded.origin_mm().map(|b| I::new(b[0], b[1]));
    let difference: [[I; 3]; 3] =
        std::array::from_fn::<_, 3, _>(|j| std::array::from_fn(|i| wanted[j][i] - actual[j][i]));
    let local = source.local_enclosure_mm().map(|b| I::new(b[0], b[1]));
    let error = std::array::from_fn(|i| {
        (0..3).fold(origin[i] - actual_origin[i], |s, j| {
            s + difference[j][i] * local[j]
        })
    });
    let position = l1(error).hi;
    let mut derivative = 0f64;
    for (t, tensor) in source.tensors().iter().enumerate() {
        for u in 0..tensor.spans()[0] {
            for v in 0..tensor.spans()[1] {
                w.step(8192)?;
                let jet = source
                    .cell(t, [u, v])?
                    .first_derivative_bounds()
                    .map_err(|_| Error::Unsupported)?;
                for d in [jet.du, jet.dv] {
                    let d = d.map(|b| I::new(b[0], b[1]));
                    derivative = derivative.max(
                        l1(std::array::from_fn(|i| {
                            (0..3).fold(I::point(0.), |s, j| s + difference[j][i] * d[j])
                        }))
                        .hi,
                    );
                }
            }
        }
    }
    let original = source.linear_columns_bounds()[0].map(|b| I::new(b[0], b[1]));
    let scale = norm(original);
    if !scale.valid() || scale.lo <= 0. {
        return Err(Error::Unsupported);
    }
    // Normalize both maps by the SAME exact original STEP scale, whose
    // enclosure is shared here. Positive scalar factors do not change unit
    // normals, and avoid overflowing cofactor squared norms at extreme scales.
    let wanted = wanted.map(|c| c.map(|v| v / scale));
    let actual = actual.map(|c| c.map(|v| v / scale));
    let cw = [
        cross(wanted[1], wanted[2]),
        cross(wanted[2], wanted[0]),
        cross(wanted[0], wanted[1]),
    ];
    let ca = [
        cross(actual[1], actual[2]),
        cross(actual[2], actual[0]),
        cross(actual[0], actual[1]),
    ];
    let frobenius = (0..3)
        .fold(I::point(0.), |s, j| {
            s + (0..3).fold(I::point(0.), |s, i| s + (cw[j][i] - ca[j][i]).square())
        })
        .sqrt();
    let normal =
        (I::point(2.) * frobenius / (I::point(1.) - I::point(placement.metric_error()))).hi;
    if !position.is_finite()
        || position > options.position_error_mm
        || !derivative.is_finite()
        || derivative > options.original_derivative_error_mm
        || !normal.is_finite()
        || normal > options.normal_chord_error
    {
        return Err(Error::Unsupported);
    }
    Ok((position, derivative, normal))
}
