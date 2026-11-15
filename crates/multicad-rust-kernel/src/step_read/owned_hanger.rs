//! Independent literal trimmed-tensor STEP readback. Parsing retains the actual
//! numeric leaves and incidence; it grants no material or Native admission.
use super::part21::{self, Entity, Value};
use crate::tapered_thread::interval::I;
use std::{
    collections::{BTreeMap, BTreeSet},
    mem::size_of,
};
mod bernstein;
mod decode;
mod lifts;
mod links;
mod numbers;
pub use lifts::{EdgeLift, LiftOptions, SurfaceLifts};
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Error {
    Invalid,
    Unsupported,
    Limit,
    Cancelled,
}
#[derive(Clone, Copy)]
pub struct ReadOptions {
    pub maximum_work: usize,
    /// Retained input, parsed arenas and decoded arrays. The bounded Part21
    /// parser has independent fixed limits; this is not an allocator peak claim.
    pub retained_budget_bytes: usize,
}
impl Default for ReadOptions {
    fn default() -> Self {
        Self {
            maximum_work: 128_000_000,
            retained_budget_bytes: 64 * 1024 * 1024,
        }
    }
}
struct Budget<'a> {
    work: usize,
    remaining: usize,
    bytes: usize,
    limit: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl Budget<'_> {
    fn charge(&mut self, n: usize) -> Result<(), Error> {
        if self.interrupt.is_some_and(|f| f()) {
            return Err(Error::Cancelled);
        }
        if n > self.remaining {
            #[cfg(test)]
            eprintln!("reader work refusal used={} next={n}", self.work);
            return Err(Error::Limit);
        }
        self.remaining = self.remaining.checked_sub(n).ok_or(Error::Limit)?;
        self.work = self.work.checked_add(n).ok_or(Error::Limit)?;
        Ok(())
    }
    fn retain(&mut self, n: usize) -> Result<(), Error> {
        self.bytes = self.bytes.checked_add(n).ok_or(Error::Limit)?;
        if self.bytes > self.limit {
            #[cfg(test)]
            eprintln!(
                "reader byte refusal bytes={} limit={} work={}",
                self.bytes, self.limit, self.work
            );
            return Err(Error::Limit);
        }
        Ok(())
    }
}
#[derive(Clone, Copy, Debug)]
pub struct Scalar<'i> {
    text: &'i str,
    value: f64,
    unit: i8,
}
impl<'i> Scalar<'i> {
    fn new(text: &'i str, value: f64) -> Self {
        let mut out = Self {
            text,
            value,
            unit: -1,
        };
        out.unit = if out.zero_spelling() {
            0
        } else if out.one_spelling() {
            1
        } else {
            -1
        };
        out
    }
    fn is_zero(&self) -> bool {
        self.unit == 0
    }
    fn is_one(&self) -> bool {
        self.unit == 1
    }
    pub fn spelling(&self) -> &'i str {
        self.text
    }
    pub fn bounds(&self) -> [f64; 2] {
        [self.value.next_down(), self.value.next_up()]
    }
    fn zero_spelling(&self) -> bool {
        if self.value != 0. {
            return false;
        }
        self.text
            .split('E')
            .next()
            .unwrap_or("")
            .bytes()
            .filter(|b| !matches!(b, b'+' | b'-' | b'.'))
            .all(|b| b == b'0')
    }
    fn one_spelling(&self) -> bool {
        if self.value != 1. {
            return false;
        }
        let text = self.text.strip_prefix('+').unwrap_or(self.text);
        if text.starts_with('-') {
            return false;
        }
        let (m, e) = text.split_once('E').unwrap_or((text, "0"));
        let Ok(e) = e.parse::<i32>() else {
            return false;
        };
        let fraction = m.split_once('.').map_or(0, |(_, f)| f.len());
        let digits = m.bytes().filter(|&b| b != b'.').collect::<Vec<_>>();
        let Some(first) = digits.iter().position(|&b| b != b'0') else {
            return false;
        };
        digits[first] == b'1'
            && digits[first + 1..].iter().all(|&b| b == b'0')
            && i64::from(e) - fraction as i64 + (digits.len() - first - 1) as i64 == 0
    }
    fn interval(&self) -> I {
        if self.is_zero() {
            I::point(0.)
        } else if self.is_one() {
            I::point(1.)
        } else {
            I::new(self.value.next_down(), self.value.next_up())
        }
    }
}
#[derive(Debug)]
pub struct Point<'i> {
    entity: u32,
    coordinates: Vec<Scalar<'i>>,
}
impl<'i> Point<'i> {
    pub fn entity(&self) -> u32 {
        self.entity
    }
    pub fn coordinates(&self) -> &[Scalar<'i>] {
        &self.coordinates
    }
}
#[derive(Debug)]
pub struct Curve<'i> {
    entity: u32,
    degree: usize,
    points: Vec<u32>,
    weights: Vec<Scalar<'i>>,
    knots: Vec<Scalar<'i>>,
}
impl<'i> Curve<'i> {
    pub fn entity(&self) -> u32 {
        self.entity
    }
    pub fn degree(&self) -> usize {
        self.degree
    }
    pub fn point_entities(&self) -> &[u32] {
        &self.points
    }
    pub fn weights(&self) -> &[Scalar<'i>] {
        &self.weights
    }
    pub fn knots(&self) -> &[Scalar<'i>] {
        &self.knots
    }
}
#[derive(Debug)]
pub struct Tensor<'i> {
    entity: u32,
    degrees: [usize; 2],
    dimensions: [usize; 2],
    points: Vec<u32>,
    weights: Vec<Scalar<'i>>,
    knots_v: Vec<Scalar<'i>>,
}
impl<'i> Tensor<'i> {
    pub fn entity(&self) -> u32 {
        self.entity
    }
    pub fn degrees(&self) -> [usize; 2] {
        self.degrees
    }
    pub fn dimensions(&self) -> [usize; 2] {
        self.dimensions
    }
    pub fn point_entities(&self) -> &[u32] {
        &self.points
    }
    pub fn weights(&self) -> &[Scalar<'i>] {
        &self.weights
    }
    pub fn knots_v(&self) -> &[Scalar<'i>] {
        &self.knots_v
    }
}
#[derive(Debug)]
pub struct Plane {
    entity: u32,
    origin: u32,
    axis: u32,
    normal: u32,
    reference: u32,
}
impl Plane {
    pub fn entity(&self) -> u32 {
        self.entity
    }
    pub fn origin_entity(&self) -> u32 {
        self.origin
    }
    pub fn axis_entity(&self) -> u32 {
        self.axis
    }
    pub fn direction_entities(&self) -> [u32; 2] {
        [self.normal, self.reference]
    }
}
#[derive(Clone, Copy, Debug)]
pub enum Surface {
    Tensor(usize),
    Plane(usize),
}
#[derive(Clone, Copy, Debug)]
pub struct Coedge {
    pub edge: u32,
    pub forward: bool,
}
#[derive(Debug)]
pub struct Face {
    entity: u32,
    surface: u32,
    kind: Surface,
    same_sense: bool,
    loops: Vec<Vec<Coedge>>,
}
impl Face {
    pub fn entity(&self) -> u32 {
        self.entity
    }
    pub fn surface_entity(&self) -> u32 {
        self.surface
    }
    pub fn surface(&self) -> Surface {
        self.kind
    }
    pub fn same_sense(&self) -> bool {
        self.same_sense
    }
    pub fn loops(&self) -> &[Vec<Coedge>] {
        &self.loops
    }
}
#[derive(Clone, Copy, Debug)]
pub struct Pcurve {
    pub surface: u32,
    pub curve: u32,
}
#[derive(Debug)]
pub struct Edge {
    entity: u32,
    vertices: [u32; 2],
    curve: u32,
    pcurves: [Pcurve; 2],
    same_sense: bool,
}
impl Edge {
    pub fn entity(&self) -> u32 {
        self.entity
    }
    pub fn vertices(&self) -> [u32; 2] {
        self.vertices
    }
    pub fn curve_entity(&self) -> u32 {
        self.curve
    }
    pub fn pcurves(&self) -> &[Pcurve; 2] {
        &self.pcurves
    }
    pub fn same_sense(&self) -> bool {
        self.same_sense
    }
}
pub struct Decoded<'i> {
    input: &'i str,
    points: BTreeMap<u32, Point<'i>>,
    curves: BTreeMap<u32, Curve<'i>>,
    tensors: Vec<Tensor<'i>>,
    planes: Vec<Plane>,
    faces: Vec<Face>,
    edges: BTreeMap<u32, Edge>,
    vertices: BTreeMap<u32, u32>,
    solid: u32,
    shell: u32,
    work: usize,
    bytes: usize,
}
impl<'i> Decoded<'i> {
    pub fn literal_input(&self) -> &'i str {
        self.input
    }
    pub fn points(&self) -> &BTreeMap<u32, Point<'i>> {
        &self.points
    }
    pub fn curves(&self) -> &BTreeMap<u32, Curve<'i>> {
        &self.curves
    }
    pub fn tensors(&self) -> &[Tensor<'i>] {
        &self.tensors
    }
    pub fn planes(&self) -> &[Plane] {
        &self.planes
    }
    pub fn faces(&self) -> &[Face] {
        &self.faces
    }
    pub fn edges(&self) -> &BTreeMap<u32, Edge> {
        &self.edges
    }
    pub fn vertices(&self) -> &BTreeMap<u32, u32> {
        &self.vertices
    }
    pub fn solid_entity(&self) -> u32 {
        self.solid
    }
    pub fn shell_entity(&self) -> u32 {
        self.shell
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.bytes
    }
}
pub fn decode<'i>(
    input: &'i str,
    options: ReadOptions,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<Decoded<'i>, Error> {
    decode::read(input, options, interrupt)
}
#[cfg(test)]
mod tests;
