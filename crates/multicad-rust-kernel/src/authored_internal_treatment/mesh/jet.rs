//! Directed second-order source jets on one closed normalized rectangle.
use crate::tapered_thread::interval::I;
#[derive(Clone, Copy)]
pub(super) struct Jet {
    pub value: I,
    pub first: [I; 2],
    pub second: [I; 3],
}
impl Jet {
    pub fn constant(value: I) -> Self {
        Self {
            value,
            first: [I::point(0.); 2],
            second: [I::point(0.); 3],
        }
    }
    pub fn variable(axis: usize) -> Self {
        let mut out = Self::constant(I::new(0., 1.));
        out.first[axis] = I::point(1.);
        out
    }
    pub fn add(self, other: Self) -> Self {
        Self {
            value: self.value + other.value,
            first: std::array::from_fn(|i| self.first[i] + other.first[i]),
            second: std::array::from_fn(|i| self.second[i] + other.second[i]),
        }
    }
    pub fn sub(self, other: Self) -> Self {
        self.add(other.scale(I::point(-1.)))
    }
    pub fn scale(self, factor: I) -> Self {
        Self {
            value: self.value * factor,
            first: self.first.map(|v| v * factor),
            second: self.second.map(|v| v * factor),
        }
    }
    pub fn mul(self, other: Self) -> Self {
        let a = self;
        let b = other;
        Self {
            value: a.value * b.value,
            first: std::array::from_fn(|i| a.first[i] * b.value + a.value * b.first[i]),
            second: [
                a.second[0] * b.value
                    + I::point(2.) * a.first[0] * b.first[0]
                    + a.value * b.second[0],
                a.second[1] * b.value
                    + a.first[0] * b.first[1]
                    + a.first[1] * b.first[0]
                    + a.value * b.second[1],
                a.second[2] * b.value
                    + I::point(2.) * a.first[1] * b.first[1]
                    + a.value * b.second[2],
            ],
        }
    }
    fn composed(self, value: I, derivative: I, second: I) -> Self {
        Self {
            value,
            first: self.first.map(|v| derivative * v),
            second: [
                derivative * self.second[0] + second * self.first[0].square(),
                derivative * self.second[1] + second * self.first[0] * self.first[1],
                derivative * self.second[2] + second * self.first[1].square(),
            ],
        }
    }
    pub fn sin(self) -> Self {
        self.composed(self.value.sin(), self.value.cos(), -self.value.sin())
    }
    pub fn cos(self) -> Self {
        self.composed(self.value.cos(), -self.value.sin(), -self.value.cos())
    }
}
