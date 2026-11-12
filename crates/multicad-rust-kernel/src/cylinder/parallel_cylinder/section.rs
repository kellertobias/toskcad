use super::*;
fn pi() -> I {
    I::new(
        std::f64::consts::PI.next_down(),
        std::f64::consts::PI.next_up(),
    )
}
fn angle(y: I, x: I) -> I {
    if x.lo > 0. {
        (y / x).atan()
    } else if x.hi < 0. {
        pi() - (y / (-x)).atan()
    } else {
        pi() * I::point(0.5) - (x / y).atan()
    }
}
fn segment(theta: I) -> Result<I, Error> {
    positive(theta)?;
    let value = if theta.hi <= 0.125 {
        let q = theta.square();
        let polynomial = I::point(2.) / I::point(3.) - q * (I::point(2.) / I::point(15.))
            + q.square() * (I::point(4.) / I::point(315.))
            - q.square() * q * (I::point(2.) / I::point(2835.));
        let tail = (theta * q.square() * q.square() * q * (I::point(4.) / I::point(155925.))).hi;
        theta * q * polynomial + I::new(0., tail)
    } else {
        theta - theta.sin() * theta.cos()
    };
    positive(value)
}
pub(super) fn area(a: Rational, b: Rational, xy: [f64; 2], d: I) -> Result<I, Error> {
    let (large, small) = if a >= b { (a, b) } else { (b, a) };
    let r = scalar(large);
    let s = scalar(small);
    let differenceq = large.checked_sub(small).ok();
    let difference = absolute(r - s);
    let contains = order(
        d - difference,
        differenceq.and_then(|q| radial_exact(xy, q)),
    )?;
    if contains != Ordering::Greater {
        return positive(pi() * s.square());
    }
    let sum = r + s;
    let external = positive(-radial_difference(
        d,
        sum,
        large.checked_add(small).ok(),
        xy,
    ))?;
    let internal = positive(radial_difference(d, difference, differenceq, xy))?;
    let chord =
        positive(external * (sum + d) * internal * (d + difference))?.sqrt() / (I::point(2.) * d);
    let x = (d + (r - s) * (r + s) / d) * I::point(0.5);
    let y = d - x;
    positive(r.square() * segment(angle(chord, x))? + s.square() * segment(angle(chord, y))?)
}
