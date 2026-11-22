//! Platform-independent transcendental functions.
//!
//! `f64::sin`, `atan2`, `hypot` and the other transcendental methods call the
//! operating system's maths library, whose last-bit rounding differs between
//! macOS, Linux and Windows. Geometry built on them then differs per platform:
//! subdivision decisions, work counts and exact readbacks diverge. These
//! methods use the pure-software `libm` crate instead, so the same document
//! evaluates to the same bits everywhere. Correctly rounded IEEE operations
//! (`sqrt`, `mul_add`, `powi`, arithmetic) are already portable and stay on std.
//! Clippy's `disallowed-methods` keeps the std variants out of geometry code.

/// Transcendental functions with identical results on every platform.
pub trait PortableMath: Sized {
    fn portable_sin(self) -> Self;
    fn portable_cos(self) -> Self;
    fn portable_sin_cos(self) -> (Self, Self);
    fn portable_tan(self) -> Self;
    fn portable_asin(self) -> Self;
    fn portable_acos(self) -> Self;
    fn portable_atan(self) -> Self;
    fn portable_atan2(self, other: Self) -> Self;
    fn portable_sinh(self) -> Self;
    fn portable_cosh(self) -> Self;
    fn portable_tanh(self) -> Self;
    fn portable_asinh(self) -> Self;
    fn portable_acosh(self) -> Self;
    fn portable_atanh(self) -> Self;
    fn portable_exp(self) -> Self;
    fn portable_exp_m1(self) -> Self;
    fn portable_ln(self) -> Self;
    fn portable_ln_1p(self) -> Self;
    fn portable_log2(self) -> Self;
    fn portable_log10(self) -> Self;
    fn portable_powf(self, exponent: Self) -> Self;
    fn portable_cbrt(self) -> Self;
    fn portable_hypot(self, other: Self) -> Self;
}

macro_rules! portable_math {
    ($float:ty, $sin:ident, $cos:ident, $sincos:ident, $tan:ident, $asin:ident, $acos:ident,
     $atan:ident, $atan2:ident, $sinh:ident, $cosh:ident, $tanh:ident, $asinh:ident,
     $acosh:ident, $atanh:ident, $exp:ident, $expm1:ident, $log:ident, $log1p:ident,
     $log2:ident, $log10:ident, $pow:ident, $cbrt:ident, $hypot:ident) => {
        impl PortableMath for $float {
            #[inline]
            fn portable_sin(self) -> Self {
                libm::$sin(self)
            }
            #[inline]
            fn portable_cos(self) -> Self {
                libm::$cos(self)
            }
            #[inline]
            fn portable_sin_cos(self) -> (Self, Self) {
                libm::$sincos(self)
            }
            #[inline]
            fn portable_tan(self) -> Self {
                libm::$tan(self)
            }
            #[inline]
            fn portable_asin(self) -> Self {
                libm::$asin(self)
            }
            #[inline]
            fn portable_acos(self) -> Self {
                libm::$acos(self)
            }
            #[inline]
            fn portable_atan(self) -> Self {
                libm::$atan(self)
            }
            #[inline]
            fn portable_atan2(self, other: Self) -> Self {
                libm::$atan2(self, other)
            }
            #[inline]
            fn portable_sinh(self) -> Self {
                libm::$sinh(self)
            }
            #[inline]
            fn portable_cosh(self) -> Self {
                libm::$cosh(self)
            }
            #[inline]
            fn portable_tanh(self) -> Self {
                libm::$tanh(self)
            }
            #[inline]
            fn portable_asinh(self) -> Self {
                libm::$asinh(self)
            }
            #[inline]
            fn portable_acosh(self) -> Self {
                libm::$acosh(self)
            }
            #[inline]
            fn portable_atanh(self) -> Self {
                libm::$atanh(self)
            }
            #[inline]
            fn portable_exp(self) -> Self {
                libm::$exp(self)
            }
            #[inline]
            fn portable_exp_m1(self) -> Self {
                libm::$expm1(self)
            }
            #[inline]
            fn portable_ln(self) -> Self {
                libm::$log(self)
            }
            #[inline]
            fn portable_ln_1p(self) -> Self {
                libm::$log1p(self)
            }
            #[inline]
            fn portable_log2(self) -> Self {
                libm::$log2(self)
            }
            #[inline]
            fn portable_log10(self) -> Self {
                libm::$log10(self)
            }
            #[inline]
            fn portable_powf(self, exponent: Self) -> Self {
                libm::$pow(self, exponent)
            }
            #[inline]
            fn portable_cbrt(self) -> Self {
                libm::$cbrt(self)
            }
            #[inline]
            fn portable_hypot(self, other: Self) -> Self {
                libm::$hypot(self, other)
            }
        }
    };
}

portable_math!(
    f64, sin, cos, sincos, tan, asin, acos, atan, atan2, sinh, cosh, tanh, asinh, acosh, atanh,
    exp, expm1, log, log1p, log2, log10, pow, cbrt, hypot
);
portable_math!(
    f32, sinf, cosf, sincosf, tanf, asinf, acosf, atanf, atan2f, sinhf, coshf, tanhf, asinhf,
    acoshf, atanhf, expf, expm1f, logf, log1pf, log2f, log10f, powf, cbrtf, hypotf
);

#[cfg(test)]
mod tests {
    use super::PortableMath;

    #[test]
    fn portable_functions_match_reference_bits_and_std_within_one_ulp() {
        // Reference bits are libm's own results; they must never depend on
        // the host maths library.
        assert_eq!(1.0_f64.portable_sin().to_bits(), 0x3feaed548f090cee);
        assert_eq!(1.0_f64.portable_cos().to_bits(), 0x3fe14a280fb5068c);
        assert_eq!(1.0_f64.portable_atan2(3.0).to_bits(), 0x3fd4978fa3269ee1);
        assert_eq!(3.0_f64.portable_hypot(4.0), 5.0);
        let ulp = |a: f64, b: f64| (a.to_bits() as i64 - b.to_bits() as i64).abs();
        for i in 1..2000 {
            let x = f64::from(i) * 0.0137 - 13.;
            let (s, c) = x.portable_sin_cos();
            assert_eq!(s.to_bits(), x.portable_sin().to_bits());
            assert_eq!(c.to_bits(), x.portable_cos().to_bits());
            #[allow(clippy::disallowed_methods)]
            {
                assert!(ulp(s, x.sin()) <= 1, "{x}");
                assert!(ulp(c, x.cos()) <= 1, "{x}");
                assert!(ulp(x.portable_atan2(0.7), x.atan2(0.7)) <= 1, "{x}");
                assert!(ulp(x.portable_hypot(0.7), x.hypot(0.7)) <= 1, "{x}");
            }
        }
    }
}
