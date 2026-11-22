//! Read-only evaluation of the wheel's generating surfaces for callers that
//! intersect, outline or section the exact tooth spaces.
use super::jet::Jet;
use super::wheel::{CUTTER_TIP, Piece, WORM_ROOT, Wheel};

impl Wheel {
    /// Sector-0 generating surface in module units over `(a, v)`: the point
    /// and its partial derivatives along `a` and `v`. `a` is the cutter
    /// radius on a flank and the counterclockwise profile parameter in
    /// `[0, 1]` on the other pieces; `v` is the wheel axial coordinate.
    pub fn surface_frame(&self, piece: Piece, a: f64, v: f64) -> [[f64; 2]; 3] {
        let [x, y] = self.surface_point(piece, Jet::variable(a, 0), Jet::variable(v, 1));
        [[x.value, y.value], [x.d[0], y.d[0]], [x.d[1], y.d[1]]]
    }

    /// Cutter radii between which the flanks are generated, in module units:
    /// the worm root and the cutter tip. A flank is retained from where it
    /// leaves the blank cylinder to the cutter tip.
    pub fn flank_radii(&self) -> [f64; 2] {
        [WORM_ROOT, CUTTER_TIP]
    }

    /// Blank cylinder radius in module units.
    pub fn blank_radius(&self) -> f64 {
        self.blank()
    }
}

#[cfg(test)]
mod tests {
    use super::super::wheel::{
        PIECES,
        tests::{spec, wheel},
    };
    use super::*;

    #[test]
    fn surface_frames_are_the_section_pieces_and_their_derivatives() {
        let wheel = wheel(spec(1., 40, 1.5, 4., 0.05));
        let [root, tip] = wheel.flank_radii();
        assert_eq!([root, tip], [3.75, 6.25]);
        assert_eq!(wheel.blank_radius(), 21.);
        for v in [-0.75, -0.2, 0.4] {
            let section = wheel.section(v);
            for piece in PIECES {
                for u in [0., 0.3, 1.] {
                    let a = match piece {
                        Piece::FlankMinus => tip - (tip - section.crossing[0]) * (1. - u),
                        Piece::FlankPlus => tip - (tip - section.crossing[1]) * u,
                        _ => u,
                    };
                    let [p, da, dv] = wheel.surface_frame(piece, a, v);
                    let exact = wheel.piece_point(&section, piece, u);
                    assert!((p[0] - exact[0]).hypot(p[1] - exact[1]) < 1e-12);
                    let h = 1e-6;
                    let [ahead, behind] =
                        [h, -h].map(|h| wheel.surface_frame(piece, a + h * (0.5 - u), v)[0]);
                    let [above, below] = [h, -h].map(|h| wheel.surface_frame(piece, a, v + h)[0]);
                    for k in 0..2 {
                        let slope = (ahead[k] - behind[k]) / (2. * h * (0.5 - u));
                        assert!((slope - da[k]).abs() < 1e-6, "{piece:?} {slope} {da:?}");
                        let rise = (above[k] - below[k]) / (2. * h);
                        assert!((rise - dv[k]).abs() < 1e-6, "{piece:?} {rise} {dv:?}");
                    }
                }
            }
        }
    }
}
