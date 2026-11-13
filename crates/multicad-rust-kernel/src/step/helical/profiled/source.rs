//! Borrowed typed geometry dispatch; no substitute recess or rebuilt thread.
use super::*;
use crate::{head_profile::Profile, phillips_drive::PhillipsScrew};
#[derive(Clone, Copy)]
pub(super) enum Kind<'a> {
    Existing(&'a Screw),
    Phillips(&'a PhillipsScrew),
}
#[derive(Clone, Copy)]
pub(crate) struct HeadSource<'a> {
    pub(super) data: outer::Data<'a>,
    pub(super) kind: Kind<'a>,
}
impl<'a> From<&'a Screw> for HeadSource<'a> {
    fn from(s: &'a Screw) -> Self {
        Self {
            data: outer::Data::from(s),
            kind: Kind::Existing(s),
        }
    }
}
impl<'a> HeadSource<'a> {
    pub fn phillips(s: &'a PhillipsScrew) -> Self {
        Self {
            data: outer::Data::from(s),
            kind: Kind::Phillips(s),
        }
    }
    pub fn existing(self) -> Option<&'a Screw> {
        match self.kind {
            Kind::Existing(s) => Some(s),
            Kind::Phillips(_) => None,
        }
    }
    pub fn phillips_geometry(self) -> Option<&'a PhillipsScrew> {
        match self.kind {
            Kind::Phillips(s) => Some(s),
            Kind::Existing(_) => None,
        }
    }
    pub fn plain_length_mm(self) -> f64 {
        self.data.plain
    }
    pub fn profile(self) -> &'a Profile {
        self.data.profile
    }
    pub fn thread(self) -> &'a Thread {
        self.data.thread
    }
    pub fn head_top_z(self) -> f64 {
        self.data.top_z()
    }
    pub fn visible_thread_interval_mm(self) -> [f64; 2] {
        self.data.visible
    }
    pub fn face_count(self) -> usize {
        match self.kind {
            Kind::Existing(s) => s.faces().len(),
            Kind::Phillips(s) => s.faces().len(),
        }
    }
}
