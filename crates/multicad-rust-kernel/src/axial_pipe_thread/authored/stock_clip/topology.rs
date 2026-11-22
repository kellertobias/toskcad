//! Physical source incidence. Three contact chart seams are intentionally
//! absent; the natural-crest/stock contact remains a real material crease.
use super::*;
#[derive(Clone, Copy, Debug)]
pub struct StockClipFaceLoop {
    face: StockClipFace,
    coedges: [StockClipCoedge; 6],
    length: usize,
}
impl StockClipFaceLoop {
    pub fn face(&self) -> StockClipFace {
        self.face
    }
    pub fn coedges(&self) -> &[StockClipCoedge] {
        &self.coedges[..self.length]
    }
}
impl StockClipPlan {
    pub fn physical_faces(&self) -> [StockClipFace; 7] {
        [
            StockClipFace::LowerFlank,
            StockClipFace::NaturalCrest,
            StockClipFace::StockCylinder,
            StockClipFace::UpperFlank,
            StockClipFace::Root,
            StockClipFace::StartCap,
            StockClipFace::EndCap,
        ]
    }
    pub fn physical_curves(&self) -> [StockClipCurve; 17] {
        // Eight split axial edges, two four-edge true end rims, and ONE
        // crease ring. All keys are the original source graph keys.
        std::array::from_fn(|i| {
            if i < 8 {
                self.curves()[i]
            } else if i < 12 {
                StockClipCurve::Ring {
                    station: StockClipStation::Start,
                    band: BANDS[i - 8],
                }
            } else if i == 12 {
                StockClipCurve::Ring {
                    station: StockClipStation::Contact,
                    band: StockClipBand::Crest,
                }
            } else {
                StockClipCurve::Ring {
                    station: StockClipStation::End,
                    band: BANDS[i - 13],
                }
            }
        })
    }
    /// One true boundary wire, with source chart seams cancelled by their
    /// identical curve owner and opposite coedge uses, never point fitting.
    pub fn face_loop(&self, face: StockClipFace) -> StockClipFaceLoop {
        let band = match face {
            StockClipFace::LowerFlank => StockClipBand::Lower,
            StockClipFace::UpperFlank => StockClipBand::Upper,
            StockClipFace::Root => StockClipBand::Root,
            _ => StockClipBand::Crest,
        };
        let before = StockClipPatch {
            slab: StockClipSlab::Before,
            band,
        }
        .coedges();
        let after = StockClipPatch {
            slab: StockClipSlab::After,
            band,
        }
        .coedges();
        let (coedges, length) = match face {
            StockClipFace::StartCap | StockClipFace::EndCap => {
                let station = if face == StockClipFace::StartCap {
                    StockClipStation::Start
                } else {
                    StockClipStation::End
                };
                // This typed branch excludes the only invalid cap station.
                let c = self
                    .cap_coedges(station)
                    .expect("true external cap station");
                ([c[0], c[1], c[2], c[3], c[0], c[0]], 4)
            }
            StockClipFace::NaturalCrest => (
                [
                    before[0], before[1], before[2], before[3], before[0], before[0],
                ],
                4,
            ),
            StockClipFace::StockCylinder => (
                [after[0], after[1], after[2], after[3], after[0], after[0]],
                4,
            ),
            _ => (
                [
                    before[0], after[0], after[1], after[2], before[2], before[3],
                ],
                6,
            ),
        };
        StockClipFaceLoop {
            face,
            coedges,
            length,
        }
    }
}
#[cfg(test)]
mod tests;
