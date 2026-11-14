//! Original numeric leaves remain authoritative; binary64 is an outward cache.
use super::*;
#[derive(Clone, Copy, Debug)]
pub struct Scalar<'input> {
    value: f64,
    spelling: &'input str,
}
impl<'input> Scalar<'input> {
    pub fn value(&self) -> f64 {
        self.value
    }
    pub fn spelling(&self) -> &'input str {
        self.spelling
    }
    /// Cost of the bounds cache query, without charging a second decimal parse.
    /// Only zero/integer authority scans original text; other finite parsed
    /// values use two constant-time adjacent-float operations.
    pub(super) fn bounds_work(&self) -> usize {
        if self.value == 0. {
            self.spelling.len() + 1
        } else if self.value.fract() == 0. && self.value.abs() <= 9_007_199_254_740_992. {
            self.spelling.len() + 32
        } else {
            1
        }
    }
    pub fn bounds(&self) -> [f64; 2] {
        // This proves an exact literal integer, rather than assuming every parsed
        // floating integer had that exact source value (e.g. 11.0000000000000001).
        if self.is_zero()
            || (self.value.fract() == 0.
                && self.value.abs() <= 9_007_199_254_740_992.
                && self.spelling == format!("{:.16E}", self.value))
        {
            [self.value, self.value]
        } else {
            [self.value.next_down(), self.value.next_up()]
        }
    }
    pub(super) fn is_zero(&self) -> bool {
        if self.value != 0. {
            return false;
        }
        self.spelling
            .split('E')
            .next()
            .unwrap()
            .bytes()
            .filter(|b| !matches!(b, b'+' | b'-' | b'.'))
            .all(|b| b == b'0')
    }
    pub(super) fn is_one(&self) -> bool {
        if self.value != 1. {
            return false;
        }
        let text = self.spelling.strip_prefix('+').unwrap_or(self.spelling);
        let (mantissa, exponent) = text.split_once('E').unwrap_or((text, "0"));
        // A <=64KiB finite decimal rounding to one has exponent magnitude
        // <=65536 (its mantissa cannot offset a larger order of magnitude).
        let exponent = exponent
            .parse::<i32>()
            .expect("unit literal exponent bounded by token length");
        let fraction = mantissa.split_once('.').map_or(0, |(_, s)| s.len());
        let mut seen = false;
        let mut trailing = 0usize;
        for digit in mantissa.bytes().filter(|&b| b != b'.') {
            if seen {
                if digit != b'0' {
                    return false;
                }
                trailing += 1;
            } else if digit != b'0' {
                if digit != b'1' {
                    return false;
                }
                seen = true;
            }
        }
        seen && i64::from(exponent) - fraction as i64 + trailing as i64 == 0
    }
}
pub(super) struct Numeric<'input> {
    tokens: Vec<&'input str>,
    next: usize,
}
impl<'input> Numeric<'input> {
    pub fn new(entity: &Entity<'input>, w: &mut Ledger<'_>) -> Result<Self> {
        // Only the actual Part21 arena reaches this module; its parser retains
        // the original source slice for every simple/complex entity component.
        let source = entity.source.expect("Part21 original source is retained");
        w.step(source.len())?;
        let tokens = part21::numeric_spellings(entity).ok_or(Error::Invalid)?;
        Ok(Self { tokens, next: 0 })
    }
    pub fn scalar(&mut self, value: &Value<'_>) -> Result<Scalar<'input>> {
        let Value::Number(value) = value else {
            return Err(Error::Invalid);
        };
        // Part21::Parser rejects nonfinite numerical values before publication.
        let spelling = *self.tokens.get(self.next).ok_or(Error::Invalid)?;
        self.next += 1;
        Ok(Scalar {
            value: *value,
            spelling,
        })
    }
    pub fn integer(&mut self, value: &Value<'_>) -> Result<u32> {
        let s = self.scalar(value)?;
        let text = s.spelling.strip_prefix('+').unwrap_or(s.spelling);
        if text.is_empty() || !text.bytes().all(|b| b.is_ascii_digit()) {
            return Err(Error::Invalid);
        }
        let text = text.trim_start_matches('0');
        if text.is_empty() {
            return Ok(0);
        }
        text.parse().map_err(|_| Error::Limit)
    }
    pub fn finish(self) -> Result<()> {
        if self.next == self.tokens.len() {
            Ok(())
        } else {
            Err(Error::Invalid)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn scalar(text: &str) -> Scalar<'_> {
        Scalar {
            value: text.parse().unwrap(),
            spelling: text,
        }
    }
    #[test]
    fn exact_unit_and_zero_tests_do_not_confuse_neighboring_decimal_values() {
        for text in ["1", "+0001.000E0", "10.E-1", "1.0000000000000000E0"] {
            assert!(scalar(text).is_one());
        }
        for text in ["1.0000000000000001", "0.99999999999999999"] {
            let s = scalar(text);
            assert_eq!(s.value(), 1.);
            assert!(!s.is_one());
            assert!(s.bounds()[0] < 1. && s.bounds()[1] > 1.);
        }
        assert!(scalar("-0.E8").is_zero());
        let tiny = scalar("1.E-1000");
        assert_eq!(tiny.value(), 0.);
        assert!(!tiny.is_zero());
        assert_ne!(tiny.bounds(), [0., 0.]);
    }
    #[test]
    fn actual_original_numeric_cursor_checks_types_consumption_and_exact_schema_integer() {
        let input = "ISO-10303-21;HEADER;FILE_DESCRIPTION(('test'),'2;1');FILE_NAME('test','',(''),(''),'test','test','');FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));ENDSEC;DATA;#1=CARTESIAN_POINT('',(0.,1.,2.));ENDSEC;END-ISO-10303-21;";
        let entities = part21::parse(input).unwrap();
        let e = &entities[&1][0];
        let mut w = Ledger {
            used: 0,
            remaining: 16_000_000,
            interrupt: None,
        };
        let n = Numeric::new(e, &mut w).unwrap();
        assert!(matches!(n.finish(), Err(Error::Invalid)));
        let mut n = Numeric::new(e, &mut w).unwrap();
        assert!(matches!(n.scalar(&Value::Enum("T")), Err(Error::Invalid)));
        let input = "ISO-10303-21;HEADER;FILE_DESCRIPTION(('test'),'2;1');FILE_NAME('test','',(''),(''),'test','test','');FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));ENDSEC;DATA;#1=B_SPLINE_CURVE_WITH_KNOTS('',0,(),.UNSPECIFIED.,.F.,.F.,(),(),.UNSPECIFIED.);ENDSEC;END-ISO-10303-21;";
        let entities = part21::parse(input).unwrap();
        let e = &entities[&1][0];
        let mut n = Numeric::new(e, &mut w).unwrap();
        assert_eq!(n.integer(&e.args[1]).unwrap(), 0);
        n.finish().unwrap();
    }
}
