//! Making floating-point values behave: traits, ordering, equality, hashing,
//! and constraints for floating-point types.
//!
//! Decorum provides proxy (wrapper) types and functions that canonicalize
//! floating-point values and provide a total ordering. This allows
//! floating-point values to be hashed and compared. Proxy types also provide
//! contraints on the values that may be represented and will panic if those
//! constraints are violated. See the [README](https://docs.rs/crate/decorum/).
//!
//! # Ordering
//!
//! `NaN` and zero are canonicalized to a single representation (called `CNaN`
//! and `C0` respectively) to provide the following total ordering for all
//! proxy types and ordering functions:
//!
//! `[ -INF | ... | C0 | ... | +INF | CNaN ]`
//!
//! Note that `NaN` is canonicalized to `CNaN`, which has a single
//! representation and supports the relations `CNaN = CNaN` and `CNaN > x | x ≠
//! CNaN`. `+0` and `-0` are also canonicalized to `C0`, which is equivalent to
//! `+0`.
//!
//! # Constraints
//!
//! The `NotNan` and `Finite` types wrap raw floating-point values and disallow
//! certain values like `NaN`, `INF`, and `-INF`. They will panic if an
//! operation or conversion invalidates these constraints. The `Ordered` type
//! allows any valid IEEE-754 value (there are no constraints). For most use
//! cases, either `Ordered` or `NotNan` are appropriate.

#![no_std]

#[cfg(feature = "std")]
extern crate std;

use core::num::FpCategory;
use core::ops::Neg;
#[allow(unused_imports)]
use num_traits::{real, Num, NumCast};

#[cfg(not(feature = "std"))]
use num_traits::float::FloatCore as Float;
#[cfg(feature = "std")]
use num_traits::Float;

// TODO: Support `f128`.

mod canonical;
mod constraint;
mod proxy;

use crate::constraint::{FiniteConstraint, NotNanConstraint, UnitConstraint};

pub use crate::canonical::{
    cmp_float, cmp_float_array, cmp_float_slice, eq_float, eq_float_array, eq_float_slice,
    hash_float, hash_float_array, hash_float_slice,
};
pub use crate::proxy::ConstrainedFloat;

/// An ordered and canonicalized floating-point value.
pub type Ordered<T> = ConstrainedFloat<T, UnitConstraint<T>>;

/// An ordered and canonicalized floating-point value that cannot be `NaN`.
///
/// If any operation results in a `NaN` value, then a panic will occur.
pub type NotNan<T> = ConstrainedFloat<T, NotNanConstraint<T>>;

/// An alias for a floating-point value that cannot be `NaN`.
pub type N32 = NotNan<f32>;
/// An alias for a floating-point value that cannot be `NaN`.
pub type N64 = NotNan<f64>;

/// An ordered and canonicalized floating-point value that must represent a
/// real number.
///
/// `NaN`, `INF`, etc. are not allowed and a panic will occur if any operation
/// results in such a value. This is sometimes referred to simply as a "real" as
/// seen in the `R32` and `R64` aliases.
pub type Finite<T> = ConstrainedFloat<T, FiniteConstraint<T>>;

/// An alias for a floating-point value that represents a real number.
///
/// The prefix "R" for "real" is used instead of "F" for "finite", because if
/// "F" were used, then this name would be very similar to `f32`,
/// differentiated only by capitalization.
pub type R32 = Finite<f32>;
/// An alias for a floating-point value that represents a real number.
///
/// The prefix "R" for "real" is used instead of "F" for "finite", because if
/// "F" were used, then this name would be very similar to `f64`,
/// differentiated only by capitalization.
pub type R64 = Finite<f64>;

/// A primitive floating-point value.
///
/// This trait differentiates types that implement floating-point traits but
/// may not be primitive types.
pub trait Primitive: Copy + Sized {}

impl Primitive for f32 {}
impl Primitive for f64 {}

/// A floating-point value that can be infinite (`-INF` or `INF`).
pub trait Infinite: Copy + NumCast {
    fn infinity() -> Self;
    fn neg_infinity() -> Self;
    fn is_infinite(self) -> bool;
    fn is_finite(self) -> bool;
}

/// A floating-point value that can be `NaN`.
pub trait Nan: Copy + NumCast {
    fn nan() -> Self;
    fn is_nan(self) -> bool;
}

/// Floating-point encoding.
///
/// Provides values and operations that directly relate to the encoding of an
/// IEEE-754 floating-point value with the exception of `-INF`, `INF`, and
/// `NaN`. See the `Infinite` and `Nan` traits.
pub trait Encoding: Copy + NumCast {
    fn max_value() -> Self;
    fn min_value() -> Self;
    fn min_positive_value() -> Self;
    fn epsilon() -> Self;
    fn classify(self) -> FpCategory;
    fn is_normal(self) -> bool;
    fn integer_decode(self) -> (u64, i16, i8);
}

/// A value that can represent a real number.
///
/// Provides values and operations that generally apply to real numbers.
pub trait Real: Copy + Neg<Output = Self> + Num + NumCast + PartialOrd {
    fn min(self, other: Self) -> Self;
    fn max(self, other: Self) -> Self;

    fn is_sign_positive(self) -> bool;
    fn is_sign_negative(self) -> bool;
    fn signum(self) -> Self;
    fn abs(self) -> Self;

    fn floor(self) -> Self;
    fn ceil(self) -> Self;
    fn round(self) -> Self;
    fn trunc(self) -> Self;
    fn fract(self) -> Self;
    fn recip(self) -> Self;

    #[cfg(feature = "std")]
    fn mul_add(self, a: Self, b: Self) -> Self;
    #[cfg(feature = "std")]
    fn abs_sub(self, other: Self) -> Self;

    #[cfg(feature = "std")]
    fn powi(self, n: i32) -> Self;
    #[cfg(feature = "std")]
    fn powf(self, n: Self) -> Self;
    #[cfg(feature = "std")]
    fn sqrt(self) -> Self;
    #[cfg(feature = "std")]
    fn cbrt(self) -> Self;
    #[cfg(feature = "std")]
    fn exp(self) -> Self;
    #[cfg(feature = "std")]
    fn exp2(self) -> Self;
    #[cfg(feature = "std")]
    fn exp_m1(self) -> Self;
    #[cfg(feature = "std")]
    fn log(self, base: Self) -> Self;
    #[cfg(feature = "std")]
    fn ln(self) -> Self;
    #[cfg(feature = "std")]
    fn log2(self) -> Self;
    #[cfg(feature = "std")]
    fn log10(self) -> Self;
    #[cfg(feature = "std")]
    fn ln_1p(self) -> Self;

    #[cfg(feature = "std")]
    fn hypot(self, other: Self) -> Self;
    #[cfg(feature = "std")]
    fn sin(self) -> Self;
    #[cfg(feature = "std")]
    fn cos(self) -> Self;
    #[cfg(feature = "std")]
    fn tan(self) -> Self;
    #[cfg(feature = "std")]
    fn asin(self) -> Self;
    #[cfg(feature = "std")]
    fn acos(self) -> Self;
    #[cfg(feature = "std")]
    fn atan(self) -> Self;
    #[cfg(feature = "std")]
    fn atan2(self, other: Self) -> Self;
    #[cfg(feature = "std")]
    fn sin_cos(self) -> (Self, Self);
    #[cfg(feature = "std")]
    fn sinh(self) -> Self;
    #[cfg(feature = "std")]
    fn cosh(self) -> Self;
    #[cfg(feature = "std")]
    fn tanh(self) -> Self;
    #[cfg(feature = "std")]
    fn asinh(self) -> Self;
    #[cfg(feature = "std")]
    fn acosh(self) -> Self;
    #[cfg(feature = "std")]
    fn atanh(self) -> Self;
}

impl<T> Infinite for T
where
    T: Float + Primitive,
{
    fn infinity() -> Self {
        Float::infinity()
    }

    fn neg_infinity() -> Self {
        Float::neg_infinity()
    }

    fn is_infinite(self) -> bool {
        Float::is_infinite(self)
    }

    fn is_finite(self) -> bool {
        Float::is_finite(self)
    }
}

impl<T> Nan for T
where
    T: Float + Primitive,
{
    fn nan() -> Self {
        Float::nan()
    }

    fn is_nan(self) -> bool {
        Float::is_nan(self)
    }
}

impl<T> Encoding for T
where
    T: Float + Primitive,
{
    fn max_value() -> Self {
        Float::max_value()
    }

    fn min_value() -> Self {
        Float::min_value()
    }

    fn min_positive_value() -> Self {
        Float::min_positive_value()
    }

    fn epsilon() -> Self {
        Float::epsilon()
    }

    fn classify(self) -> FpCategory {
        Float::classify(self)
    }

    fn is_normal(self) -> bool {
        Float::is_normal(self)
    }

    fn integer_decode(self) -> (u64, i16, i8) {
        Float::integer_decode(self)
    }
}

impl<T> Real for T
where
    T: Float + Primitive,
{
    fn min(self, other: Self) -> Self {
        Float::min(self, other)
    }

    fn max(self, other: Self) -> Self {
        Float::max(self, other)
    }

    fn is_sign_positive(self) -> bool {
        Float::is_sign_positive(self)
    }

    fn is_sign_negative(self) -> bool {
        Float::is_sign_negative(self)
    }

    fn signum(self) -> Self {
        Float::signum(self)
    }

    fn abs(self) -> Self {
        Float::abs(self)
    }

    fn floor(self) -> Self {
        Float::floor(self)
    }

    fn ceil(self) -> Self {
        Float::ceil(self)
    }

    fn round(self) -> Self {
        Float::round(self)
    }

    fn trunc(self) -> Self {
        Float::trunc(self)
    }

    fn fract(self) -> Self {
        Float::fract(self)
    }

    fn recip(self) -> Self {
        Float::recip(self)
    }

    #[cfg(feature = "std")]
    fn mul_add(self, a: Self, b: Self) -> Self {
        Float::mul_add(self, a, b)
    }

    #[cfg(feature = "std")]
    fn abs_sub(self, other: Self) -> Self {
        Float::abs_sub(self, other)
    }

    #[cfg(feature = "std")]
    fn powi(self, n: i32) -> Self {
        Float::powi(self, n)
    }

    #[cfg(feature = "std")]
    fn powf(self, n: Self) -> Self {
        Float::powf(self, n)
    }

    #[cfg(feature = "std")]
    fn sqrt(self) -> Self {
        Float::sqrt(self)
    }

    #[cfg(feature = "std")]
    fn cbrt(self) -> Self {
        Float::cbrt(self)
    }

    #[cfg(feature = "std")]
    fn exp(self) -> Self {
        Float::exp(self)
    }

    #[cfg(feature = "std")]
    fn exp2(self) -> Self {
        Float::exp2(self)
    }

    #[cfg(feature = "std")]
    fn exp_m1(self) -> Self {
        Float::exp_m1(self)
    }

    #[cfg(feature = "std")]
    fn log(self, base: Self) -> Self {
        Float::log(self, base)
    }

    #[cfg(feature = "std")]
    fn ln(self) -> Self {
        Float::ln(self)
    }

    #[cfg(feature = "std")]
    fn log2(self) -> Self {
        Float::log2(self)
    }

    #[cfg(feature = "std")]
    fn log10(self) -> Self {
        Float::log10(self)
    }

    #[cfg(feature = "std")]
    fn ln_1p(self) -> Self {
        Float::ln_1p(self)
    }

    #[cfg(feature = "std")]
    fn hypot(self, other: Self) -> Self {
        Float::hypot(self, other)
    }

    #[cfg(feature = "std")]
    fn sin(self) -> Self {
        Float::sin(self)
    }

    #[cfg(feature = "std")]
    fn cos(self) -> Self {
        Float::cos(self)
    }

    #[cfg(feature = "std")]
    fn tan(self) -> Self {
        Float::tan(self)
    }

    #[cfg(feature = "std")]
    fn asin(self) -> Self {
        Float::asin(self)
    }

    #[cfg(feature = "std")]
    fn acos(self) -> Self {
        Float::acos(self)
    }

    #[cfg(feature = "std")]
    fn atan(self) -> Self {
        Float::atan(self)
    }

    #[cfg(feature = "std")]
    fn atan2(self, other: Self) -> Self {
        Float::atan2(self, other)
    }

    #[cfg(feature = "std")]
    fn sin_cos(self) -> (Self, Self) {
        Float::sin_cos(self)
    }

    #[cfg(feature = "std")]
    fn sinh(self) -> Self {
        Float::sinh(self)
    }

    #[cfg(feature = "std")]
    fn cosh(self) -> Self {
        Float::cosh(self)
    }

    #[cfg(feature = "std")]
    fn tanh(self) -> Self {
        Float::tanh(self)
    }

    #[cfg(feature = "std")]
    fn asinh(self) -> Self {
        Float::asinh(self)
    }

    #[cfg(feature = "std")]
    fn acosh(self) -> Self {
        Float::acosh(self)
    }

    #[cfg(feature = "std")]
    fn atanh(self) -> Self {
        Float::atanh(self)
    }
}

/// Implements the `Real` trait from
/// [num-traits](https://crates.io/crates/num-traits) in terms of Decorum's
/// numeric traits. Does nothing if the `std` feature is disabled.
///
/// This is not generic, because the blanket implementation provided by
/// num-traits prevents a constraint-based implementation. Instead, this macro
/// must be applied manually to each proxy type exported by Decorum that is
/// `Real` but not `Float`.
///
/// See the following issues:
///
///   https://github.com/olson-sean-k/decorum/issues/10
///   https://github.com/rust-num/num-traits/issues/49
macro_rules! real {
    (proxy => $T:ty) => {
        #[cfg(feature = "std")]
        impl real::Real for $T {
            fn max_value() -> Self {
                Encoding::max_value()
            }

            fn min_value() -> Self {
                Encoding::min_value()
            }

            fn min_positive_value() -> Self {
                Encoding::min_positive_value()
            }

            fn epsilon() -> Self {
                Encoding::epsilon()
            }

            fn min(self, other: Self) -> Self {
                Real::min(self, other)
            }

            fn max(self, other: Self) -> Self {
                Real::max(self, other)
            }

            fn is_sign_positive(self) -> bool {
                Real::is_sign_positive(self)
            }

            fn is_sign_negative(self) -> bool {
                Real::is_sign_negative(self)
            }

            fn signum(self) -> Self {
                Real::signum(self)
            }

            fn abs(self) -> Self {
                Real::abs(self)
            }

            fn floor(self) -> Self {
                Real::floor(self)
            }

            fn ceil(self) -> Self {
                Real::ceil(self)
            }

            fn round(self) -> Self {
                Real::round(self)
            }

            fn trunc(self) -> Self {
                Real::trunc(self)
            }

            fn fract(self) -> Self {
                Real::fract(self)
            }

            fn recip(self) -> Self {
                Real::recip(self)
            }

            fn mul_add(self, a: Self, b: Self) -> Self {
                Real::mul_add(self, a, b)
            }

            fn abs_sub(self, other: Self) -> Self {
                Real::abs_sub(self, other)
            }

            fn powi(self, n: i32) -> Self {
                Real::powi(self, n)
            }

            fn powf(self, n: Self) -> Self {
                Real::powf(self, n)
            }

            fn sqrt(self) -> Self {
                Real::sqrt(self)
            }

            fn cbrt(self) -> Self {
                Real::cbrt(self)
            }

            fn exp(self) -> Self {
                Real::exp(self)
            }

            fn exp2(self) -> Self {
                Real::exp2(self)
            }

            fn exp_m1(self) -> Self {
                Real::exp_m1(self)
            }

            fn log(self, base: Self) -> Self {
                Real::log(self, base)
            }

            fn ln(self) -> Self {
                Real::ln(self)
            }

            fn log2(self) -> Self {
                Real::log2(self)
            }

            fn log10(self) -> Self {
                Real::log10(self)
            }

            fn to_degrees(self) -> Self {
                Self::from_inner(self.into_inner().to_degrees())
            }

            fn to_radians(self) -> Self {
                Self::from_inner(self.into_inner().to_radians())
            }

            fn ln_1p(self) -> Self {
                Real::ln_1p(self)
            }

            fn hypot(self, other: Self) -> Self {
                Real::hypot(self, other)
            }

            fn sin(self) -> Self {
                Real::sin(self)
            }

            fn cos(self) -> Self {
                Real::cos(self)
            }

            fn tan(self) -> Self {
                Real::tan(self)
            }

            fn asin(self) -> Self {
                Real::asin(self)
            }

            fn acos(self) -> Self {
                Real::acos(self)
            }

            fn atan(self) -> Self {
                Real::atan(self)
            }

            fn atan2(self, other: Self) -> Self {
                Real::atan2(self, other)
            }

            fn sin_cos(self) -> (Self, Self) {
                Real::sin_cos(self)
            }

            fn sinh(self) -> Self {
                Real::sinh(self)
            }

            fn cosh(self) -> Self {
                Real::cosh(self)
            }

            fn tanh(self) -> Self {
                Real::tanh(self)
            }

            fn asinh(self) -> Self {
                Real::asinh(self)
            }

            fn acosh(self) -> Self {
                Real::acosh(self)
            }

            fn atanh(self) -> Self {
                Real::atanh(self)
            }
        }
    };
}
real!(proxy => N32);
real!(proxy => N64);
real!(proxy => R32);
real!(proxy => R64);
