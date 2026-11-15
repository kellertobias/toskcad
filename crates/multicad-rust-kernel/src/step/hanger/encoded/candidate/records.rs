use super::super::decimal::DecimalCache;
use super::*;
use std::rc::Rc;
pub(super) struct Records<'a> {
    file: StepFile,
    pub(super) w: Work<'a>,
    m: Memory,
    cache: DecimalCache,
    literals: super::numbers::Literals,
    points: super::numbers::Points,
    options: FileOptions,
    interrupt: Option<&'a dyn Fn() -> bool>,
    stages: usize,
    total: usize,
    maximum: usize,
    peak: usize,
}
struct Count(usize);
impl Write for Count {
    fn write_str(&mut self, s: &str) -> fmt::Result {
        self.0 = self.0.checked_add(s.len()).ok_or(fmt::Error)?;
        Ok(())
    }
}
struct Texts<'a>(&'a [Rc<FiniteDecimal>]);
impl fmt::Display for Texts<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("(")?;
        for (i, x) in self.0.iter().enumerate() {
            if i > 0 {
                f.write_str(",")?
            }
            f.write_str(x.as_str())?;
        }
        f.write_str(")")
    }
}
struct Multiplicities(usize, usize);
impl fmt::Display for Multiplicities {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("(")?;
        for i in 0..self.0 {
            if i > 0 {
                f.write_str(",")?
            }
            write!(f, "{}", self.1 + usize::from(i == 0 || i + 1 == self.0))?;
        }
        f.write_str(")")
    }
}
struct Rows<'a>(&'a [Entity], usize);
impl fmt::Display for Rows<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("(")?;
        for (i, row) in self.0.chunks(self.1).enumerate() {
            if i > 0 {
                f.write_str(",")?
            }
            write!(f, "{}", Entities(row))?;
        }
        f.write_str(")")
    }
}
struct Weights<'a>(&'a [Rc<FiniteDecimal>], usize);
impl fmt::Display for Weights<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("(")?;
        for (i, w) in self.0.iter().enumerate() {
            if i > 0 {
                f.write_str(",")?
            }
            f.write_str("(")?;
            for j in 0..self.1 {
                if j > 0 {
                    f.write_str(",")?
                }
                f.write_str(w.as_str())?;
            }
            f.write_str(")")?;
        }
        f.write_str(")")
    }
}
impl<'a> Records<'a> {
    pub(super) fn new(
        options: FileOptions,
        inherited: usize,
        interrupt: Option<&'a dyn Fn() -> bool>,
    ) -> Result<Self, Error> {
        if options.maximum_encoding_work == 0
            || options.maximum_encoding_work > 256_000_000
            || options.maximum_stage_work == 0
            || options.maximum_stage_work > 16_000_000
        {
            return Err(Error::Invalid);
        }
        let domain = DomainOptions {
            maximum_work: options
                .maximum_stage_work
                .min(options.maximum_encoding_work),
            maximum_bytes: options.maximum_bytes,
        };
        let (mut w, mut m) = crate::step::hanger::domain::ledgers(domain, inherited, interrupt)?;
        m.add(SCRATCH + size_of::<Self>())?;
        let file = StepFile::new();
        m.add(file.data.capacity())?;
        let cache = DecimalCache::new(&mut m, &mut w)?;
        let literals = super::numbers::Literals::new(&mut m)?;
        let points = super::numbers::Points::new(&mut m)?;
        let peak = m.used;
        Ok(Self {
            file,
            w,
            m,
            cache,
            literals,
            points,
            options,
            interrupt,
            stages: 0,
            total: 0,
            maximum: 0,
            peak,
        })
    }
    pub(super) fn finish_stage(&mut self) -> Result<(), Error> {
        self.total = self.total.checked_add(self.w.used).ok_or(Error::Limit)?;
        self.maximum = self.maximum.max(self.w.used);
        self.stages = self.stages.checked_add(1).ok_or(Error::Limit)?;
        let remaining = self
            .options
            .maximum_encoding_work
            .checked_sub(self.total)
            .ok_or(Error::Limit)?;
        self.w = Work::new(
            self.options.maximum_stage_work.min(remaining),
            self.interrupt,
        )?;
        self.w.charge(1)?;
        Ok(())
    }
    pub(super) fn hold_bytes(&mut self, n: usize) -> Result<(), Error> {
        self.m.add(n)?;
        self.peak = self.peak.max(self.m.used);
        Ok(())
    }
    pub(super) fn reserve_indices(&mut self, n: usize) -> Result<(), Error> {
        self.w.charge(n)?;
        self.m.add(n.checked_mul(64).ok_or(Error::Limit)?)?;
        self.peak = self.peak.max(self.m.used);
        Ok(())
    }
    fn room(&mut self, n: usize) -> Result<(), Error> {
        let required = self.file.data.len().checked_add(n).ok_or(Error::Limit)?;
        if required > self.file.data.capacity() {
            let old = self.file.data.capacity();
            let target = required.max(if old < 8 * 1024 * 1024 {
                old.checked_mul(2).ok_or(Error::Limit)?
            } else {
                old.checked_add(1024 * 1024).ok_or(Error::Limit)?
            });
            self.m.add(target - old)?;
            if self.file.data.len() > self.w.remaining() {
                self.finish_stage()?;
            }
            self.w.charge(self.file.data.len())?;
            self.file
                .data
                .try_reserve_exact(target - self.file.data.len())
                .map_err(|_| Error::Limit)?;
            self.m.add(self.file.data.capacity() - target)?;
            self.peak = self.peak.max(self.m.used);
        }
        Ok(())
    }
    pub(super) fn entity(&mut self, args: fmt::Arguments<'_>) -> Result<Entity, Error> {
        let mut count = Count(0);
        fmt::write(&mut count, args).map_err(|_| Error::Limit)?;
        // Entity index is bounded by emitted bytes, hence at most20 decimal digits.
        let bytes = count.0.checked_add(24).ok_or(Error::Limit)?;
        self.w.charge(bytes)?;
        self.room(bytes)?;
        Ok(self.file.entity(args))
    }
    fn format(&mut self, q: &N) -> Result<Rc<FiniteDecimal>, Error> {
        self.literals
            .format(q, &mut self.cache, &mut self.m, &mut self.w, &mut self.peak)
    }
    fn release(&mut self, values: Vec<Rc<FiniteDecimal>>) {
        super::numbers::release(values, &mut self.m)
    }
    pub(super) fn point(&mut self, point: &[N], dim: usize) -> Result<Entity, Error> {
        let lookup = self.points.lookup(&point[..dim], &mut self.w)?;
        if let Some(entity) = lookup.entity {
            return Ok(entity);
        }
        let values = point[..dim]
            .iter()
            .map(|q| self.format(q))
            .collect::<Result<Vec<_>, _>>()?;
        let result = self.entity(format_args!("CARTESIAN_POINT('',{})", Texts(&values)));
        self.release(values);
        let entity = result?;
        self.points
            .insert(lookup, entity, &mut self.m, &mut self.w)?;
        self.peak = self.peak.max(self.m.used);
        Ok(entity)
    }
    pub(super) fn encoded_point(
        &mut self,
        point: &[EncodedScalar; 3],
        dim: usize,
        shift: Option<&N>,
    ) -> Result<Entity, Error> {
        let mut values = point
            .iter()
            .map(|q| q.number(&mut self.w))
            .collect::<Result<Vec<_>, _>>()?;
        if let Some(shift) = shift {
            values[1] = values[1].add(shift, &mut self.w)?;
        }
        self.point(&values, dim)
    }
    pub(super) fn binary_knots(&mut self, knots: &[f64]) -> Result<Vec<N>, Error> {
        knots.iter().map(|&q| binary(q, &mut self.w)).collect()
    }
    fn knot_texts(&mut self, knots: &[N]) -> Result<Vec<Rc<FiniteDecimal>>, Error> {
        knots.iter().map(|q| self.format(q)).collect()
    }
    pub(super) fn curve(
        &mut self,
        points: &[[EncodedScalar; 3]],
        weights: &[EncodedScalar],
        knots: &[N],
        degree: usize,
        dim: usize,
        shift: Option<N>,
    ) -> Result<Entity, Error> {
        let mut poles = Vec::with_capacity(points.len());
        for (i, p) in points.iter().enumerate() {
            if i % 64 == 0 {
                self.finish_stage()?;
            }
            poles.push(self.encoded_point(p, dim, shift.as_ref())?);
        }
        let mut ws = Vec::with_capacity(points.len());
        for i in 0..points.len() {
            let q = if weights.is_empty() {
                N::integer(1)
            } else {
                weights[if weights.len() == 1 { 0 } else { i }].number(&mut self.w)?
            };
            ws.push(self.format(&q)?);
        }
        let ks = self.knot_texts(knots)?;
        let result=self.entity(format_args!("(BOUNDED_CURVE()B_SPLINE_CURVE({degree},{},.UNSPECIFIED.,.F.,.F.)B_SPLINE_CURVE_WITH_KNOTS({},{},.UNSPECIFIED.)CURVE()GEOMETRIC_REPRESENTATION_ITEM()RATIONAL_B_SPLINE_CURVE({})REPRESENTATION_ITEM(''))",Entities(&poles),Multiplicities(knots.len(),degree),Texts(&ks),Texts(&ws)));
        self.release(ks);
        self.release(ws);
        result
    }
    pub(super) fn line(&mut self, points: [[N; 2]; 2], knots: &[N; 2]) -> Result<Entity, Error> {
        let poles = [self.point(&points[0], 2)?, self.point(&points[1], 2)?];
        let ks = self.knot_texts(knots)?;
        let result = self.entity(format_args!(
            "B_SPLINE_CURVE_WITH_KNOTS('',1,{},.UNSPECIFIED.,.F.,.F.,(2,2),{},.UNSPECIFIED.)",
            Entities(&poles),
            Texts(&ks)
        ));
        self.release(ks);
        result
    }
    pub(super) fn surface(
        &mut self,
        points: &[[EncodedScalar; 3]],
        dimensions: [usize; 2],
        degrees: [usize; 2],
        weights: &[EncodedScalar],
        knots: &[f64],
    ) -> Result<Entity, Error> {
        let mut poles = Vec::with_capacity(points.len());
        for (i, p) in points.iter().enumerate() {
            if i % 64 == 0 {
                self.finish_stage()?;
            }
            poles.push(self.encoded_point(p, 3, None)?);
        }
        let mut ws = Vec::with_capacity(dimensions[0]);
        for i in 0..dimensions[0] {
            let q = if weights.is_empty() {
                N::integer(1)
            } else {
                weights[i].number(&mut self.w)?
            };
            ws.push(self.format(&q)?);
        }
        let knots = self.binary_knots(knots)?;
        let ks = self.knot_texts(&knots)?;
        let [du, dv] = degrees;
        let result=self.entity(format_args!("(BOUNDED_SURFACE()B_SPLINE_SURFACE({du},{dv},{},.UNSPECIFIED.,.F.,.F.,.F.)B_SPLINE_SURFACE_WITH_KNOTS(({},{}),{},(0.,1.),{},.UNSPECIFIED.)GEOMETRIC_REPRESENTATION_ITEM()RATIONAL_B_SPLINE_SURFACE({})REPRESENTATION_ITEM('')SURFACE())",Rows(&poles,dimensions[1]),du+1,du+1,Multiplicities(knots.len(),dv),Texts(&ks),Weights(&ws,dimensions[1])));
        self.release(ks);
        self.release(ws);
        result
    }
    pub(super) fn pcurve(
        &mut self,
        surface: Entity,
        curve: Entity,
        context: Entity,
    ) -> Result<Entity, Error> {
        let definition = self.entity(format_args!(
            "DEFINITIONAL_REPRESENTATION('',({curve}),{context})"
        ))?;
        self.entity(format_args!("PCURVE('',{surface},{definition})"))
    }
    pub(super) fn finish(mut self, solid: Entity, known: usize) -> Result<StepCandidate, Error> {
        // Fixed STEP footer from the shared writer contains only static metadata and
        // one solid reference;8KiB bounds all of its entities and temporary lists.
        self.w.charge(8192)?;
        self.room(8192)?;
        let before = self.file.data.len();
        self.file.finish(vec![solid]);
        if self.file.data.len() - before > 8192 {
            return Err(Error::Limit);
        }
        self.total = self.total.checked_add(self.w.used).ok_or(Error::Limit)?;
        self.maximum = self.maximum.max(self.w.used);
        self.stages = self.stages.checked_add(1).ok_or(Error::Limit)?;
        let known_work = known.checked_add(self.total).ok_or(Error::Limit)?;
        let bytes = size_of::<StepCandidate>() + self.file.data.capacity();
        Ok(StepCandidate {
            text: self.file.data,
            encoding_work: self.total,
            maximum_stage_work: self.maximum,
            stages: self.stages,
            known_work,
            bytes,
            peak: self.peak,
        })
    }
}
