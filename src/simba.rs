mod simba_impl {
    use simba::scalar::{RealField, SubsetOf};
    use num_traits::{Zero, One};
    use crate::NotNan;

    use crate::Real;

    impl<T: simba::simd::PrimitiveSimdValue> simba::simd::PrimitiveSimdValue for crate::NotNan<T>
    where
        T: simba::simd::SimdValue
            + crate::Primitive
            + crate::Infinite
            + crate::IntrinsicOrd
            + crate::Nan
            + crate::Real
    {}

    impl<T> simba::simd::SimdValue for crate::NotNan<T>
    where
        T: simba::simd::SimdValue
            + crate::Primitive
            + crate::Infinite
            + crate::IntrinsicOrd
            + crate::Nan
            + crate::Real
    {
        type Element = Self;
        type SimdBool = bool;

        #[inline(always)]
        fn lanes() -> usize {
            1
        }

        #[inline(always)]
        fn splat(val: Self::Element) -> Self {
            val
        }

        #[inline(always)]
        fn extract(&self, _: usize) -> Self::Element {
            *self
        }

        #[inline(always)]
        unsafe fn extract_unchecked(&self, _: usize) -> Self::Element {
            *self
        }

        #[inline(always)]
        fn replace(&mut self, _: usize, val: Self::Element) {
            *self = val
        }

        #[inline(always)]
        unsafe fn replace_unchecked(&mut self, _: usize, val: Self::Element) {
            *self = val
        }

        #[inline(always)]
        fn select(self, cond: Self::SimdBool, other: Self) -> Self {
            if cond {
                self
            } else {
                other
            }
        }
    }

    impl<T> simba::scalar::RealField for crate::NotNan<T>
    where
        f64: simba::scalar::SubsetOf<T>,
        T: simba::scalar::RealField
            + crate::Primitive
            + crate::Infinite
            + crate::IntrinsicOrd
            + crate::Nan
            + crate::Real
    {
        #[inline]
        fn is_sign_positive(&self) -> bool {
            crate::Encoding::is_sign_positive(self.clone())
        }

        #[inline]
        fn is_sign_negative(&self) -> bool {
            crate::Encoding::is_sign_negative(self.clone())
        }

        #[inline(always)]
        fn copysign(self, to: Self) -> Self {
            // let signbit = (-T::zero()).to_bits();
            // Self::from_bits((signbit & self.to_bits()) | ((!signbit) & to.to_bits()))
            todo!()
        }

        #[inline]
        fn max(self, other: Self) -> Self {
            std::cmp::Ord::max(self, other)
        }

        #[inline]
        fn min(self, other: Self) -> Self {
            std::cmp::Ord::min(self, other)
        }

        #[inline]
        fn clamp(self, min: Self, max: Self) -> Self {
            if self < min {
                min
            } else if self > max {
                max
            } else {
                self
            }
        }

        #[inline]
        fn atan2(self, other: Self) -> Self {
            crate::Real::atan2(self, other)
        }

        /// Archimedes' constant.
        #[inline]
        fn pi() -> Self {
            crate::Real::PI
        }

        /// 2.0 * pi.
        #[inline]
        fn two_pi() -> Self {
            crate::Real::TAU
        }

        /// pi / 2.0.
        #[inline]
        fn frac_pi_2() -> Self {
            crate::Real::FRAC_PI_2
        }

        /// pi / 3.0.
        #[inline]
        fn frac_pi_3() -> Self {
            crate::Real::FRAC_PI_3
        }

        /// pi / 4.0.
        #[inline]
        fn frac_pi_4() -> Self {
            crate::Real::FRAC_PI_4
        }

        /// pi / 6.0.
        #[inline]
        fn frac_pi_6() -> Self {
            crate::Real::FRAC_PI_6
        }

        /// pi / 8.0.
        #[inline]
        fn frac_pi_8() -> Self {
            crate::Real::FRAC_PI_8
        }

        /// 1.0 / pi.
        #[inline]
        fn frac_1_pi() -> Self {
            crate::Real::FRAC_1_PI
        }

        /// 2.0 / pi.
        #[inline]
        fn frac_2_pi() -> Self {
            crate::Real::FRAC_2_PI
        }

        /// 2.0 / sqrt(pi).
        #[inline]
        fn frac_2_sqrt_pi() -> Self {
            crate::Real::FRAC_2_SQRT_PI
        }


        /// Euler's number.
        #[inline]
        fn e() -> Self {
            crate::Real::E
        }

        /// log2(e).
        #[inline]
        fn log2_e() -> Self {
            crate::Real::LOG2_E
        }

        /// log10(e).
        #[inline]
        fn log10_e() -> Self {
            crate::Real::LOG10_E
        }

        /// ln(2.0).
        #[inline]
        fn ln_2() -> Self {
            crate::Real::LN_2
        }

        /// ln(10.0).
        #[inline]
        fn ln_10() -> Self {
            crate::Real::LN_10
        }

        fn min_value() -> Option<Self> {
            T::min_value().and_then(|t| NotNan::new(t).ok())
        }

        fn max_value() -> Option<Self> {
            T::max_value().and_then(|t| NotNan::new(t).ok())
        }
    }


    impl<T> simba::scalar::SubsetOf<Self> for crate::NotNan<T>
    where
        T: simba::scalar::SubsetOf<T>
            + crate::Primitive
            + crate::Infinite
            + crate::IntrinsicOrd
            + crate::Nan
            + crate::Real
    {
        #[inline]
        fn to_superset(&self) -> Self {
            *self
        }

        #[inline]
        fn from_superset_unchecked(element: &Self) -> Self {
            *element
        }

        #[inline]
        fn is_in_subset(_: &Self) -> bool {
            true
        }
    }

    impl<T> simba::scalar::SupersetOf<f64> for crate::NotNan<T>
    where
        f64: simba::scalar::SubsetOf<T>,
        T: crate::Primitive
            + crate::Infinite
            + crate::IntrinsicOrd
            + crate::Nan
            + crate::Real
    {
        #[inline]
        fn to_subset(&self) -> Option<f64> {
            f64::from_superset(&self.into_inner())
        }

        #[inline]
        fn is_in_subset(&self) -> bool {
            <f64 as simba::scalar::SubsetOf::<T>>::is_in_subset(&self.into_inner())
        }

        #[inline]
        fn to_subset_unchecked(&self) -> f64 {
            f64::from_superset_unchecked(&self.into_inner())
        }

        #[inline]
        fn from_subset(element: &f64) -> Self {
            NotNan::new(element.to_superset()).unwrap()
        }
    }

    impl<T> simba::scalar::SupersetOf<f32> for crate::NotNan<T>
    where
        f32: simba::scalar::SubsetOf<T>,
        T: crate::Primitive
            + crate::Infinite
            + crate::IntrinsicOrd
            + crate::Nan
            + crate::Real
    {
        #[inline]
        fn to_subset(&self) -> Option<f32> {
            f32::from_superset(&self.into_inner())
        }

        #[inline]
        fn is_in_subset(&self) -> bool {
            <f32 as simba::scalar::SubsetOf::<T>>::is_in_subset(&self.into_inner())
        }

        #[inline]
        fn to_subset_unchecked(&self) -> f32 {
            f32::from_superset_unchecked(&self.into_inner())
        }

        #[inline]
        fn from_subset(element: &f32) -> Self {
            NotNan::new(element.to_superset()).unwrap()
        }
    }

    impl<T> simba::scalar::SubsetOf<T> for crate::NotNan<T>
    where T: crate::Primitive
            + crate::Infinite
            + crate::IntrinsicOrd
            + crate::Nan
            + crate::Real
    {
        #[inline]
        fn to_superset(&self) -> T {
            self.into_inner()
        }

        #[inline]
        fn from_superset_unchecked(element: &T) -> Self {
            Self::new_unchecked(*element)
        }

        #[inline]
        fn is_in_subset(element: &T) -> bool {
            !element.is_nan()
        }
    }

    impl<T> simba::scalar::Field for crate::NotNan<T>
    where
        T: simba::scalar::Field
            + crate::Primitive
            + crate::Infinite
            + crate::IntrinsicOrd
            + crate::Nan
            + crate::Real
    {}

    impl<T> simba::scalar::ComplexField for crate::NotNan<T>
    where
        f64: simba::scalar::SubsetOf<T>,
        T: simba::scalar::RealField
            + simba::scalar::ComplexField
            + crate::Primitive
            + crate::Infinite
            + crate::IntrinsicOrd
            + crate::Nan
            + crate::Real
    {
        type RealField = Self;

        #[inline]
        fn from_real(re: Self::RealField) -> Self {
            re
        }

        #[inline]
        fn real(self) -> Self::RealField {
            self
        }

        #[inline]
        fn imaginary(self) -> Self::RealField {
            Self::zero()
        }

        #[inline]
        fn norm1(self) -> Self::RealField {
            Self::abs(self)
        }

        #[inline]
        fn modulus(self) -> Self::RealField {
            Self::abs(self)
        }

        #[inline]
        fn modulus_squared(self) -> Self::RealField {
            self * self
        }

        #[inline]
        fn argument(self) -> Self::RealField {
            if self >= Self::zero() {
                Self::zero()
            } else {
                Self::pi()
            }
        }

        #[inline]
        fn to_exp(self) -> (Self, Self) {
            if self >= Self::zero() {
                (self, Self::one())
            } else {
                (-self, -Self::one())
            }
        }

        #[inline]
        fn recip(self) -> Self {
            crate::Real::recip(self)
        }

        #[inline]
        fn conjugate(self) -> Self {
            self
        }

        #[inline]
        fn scale(self, factor: Self::RealField) -> Self {
            self * factor
        }

        #[inline]
        fn unscale(self, factor: Self::RealField) -> Self {
            self / factor
        }

        #[inline]
        fn floor(self) -> Self {
            crate::Real::floor(self)
        }

        #[inline]
        fn ceil(self) -> Self {
            crate::Real::ceil(self)
        }

        #[inline]
        fn round(self) -> Self {
            crate::Real::round(self)
        }

        #[inline]
        fn trunc(self) -> Self {
            crate::Real::trunc(self)
        }

        #[inline]
        fn fract(self) -> Self {
            crate::Real::fract(self)
        }

        #[inline]
        fn abs(self) -> Self {
            crate::Signed::abs(&self)
        }

        #[inline]
        fn signum(self) -> Self {
            crate::Signed::signum(&self)
        }

        #[inline]
        fn mul_add(self, a: Self, b: Self) -> Self {
            crate::Real::mul_add(self, a, b)
        }

        #[cfg(feature = "std")]
        #[inline]
        fn powi(self, n: i32) -> Self {
            crate::Real::powi(self, n)
        }

        #[cfg(not(feature = "std"))]
        #[inline]
        fn powi(self, n: i32) -> Self {
            // FIXME: is there a more accurate solution?
            Self::powf(self, n as Self)
        }

        #[inline]
        fn powf(self, n: Self) -> Self {
            Real::powf(self, n)
        }

        #[inline]
        fn powc(self, n: Self) -> Self {
            // Same as powf.
            Real::powf(self, n)
        }

        #[inline]
        fn sqrt(self) -> Self {
            Real::sqrt(self)
        }

        #[inline]
        fn try_sqrt(self) -> Option<Self> {
            if self >= Self::zero() {
                Some(Real::sqrt(self))
            } else {
                None
            }
        }

        #[inline]
        fn exp(self) -> Self {
            Real::exp(self)
        }

        #[inline]
        fn exp2(self) -> Self {
            Real::exp2(self)
        }


        #[inline]
        fn exp_m1(self) -> Self {
            Real::exp_m1(self)
        }

        #[inline]
        fn ln_1p(self) -> Self {
            Real::ln_1p(self)
        }

        #[inline]
        fn ln(self) -> Self {
            Real::ln(self)
        }

        #[inline]
        fn log(self, base: Self) -> Self {
            Real::log(self, base)
        }

        #[inline]
        fn log2(self) -> Self {
            Real::log2(self)
        }

        #[inline]
        fn log10(self) -> Self {
            Real::log10(self)
        }

        #[inline]
        fn cbrt(self) -> Self {
            Real::cbrt(self)
        }

        #[inline]
        fn hypot(self, other: Self) -> Self::RealField {
            Real::hypot(self, other)
        }

        #[inline]
        fn sin(self) -> Self {
            Real::sin(self)
        }

        #[inline]
        fn cos(self) -> Self {
            Real::cos(self)
        }

        #[inline]
        fn tan(self) -> Self {
            Real::tan(self)
        }

        #[inline]
        fn asin(self) -> Self {
            Real::asin(self)
        }

        #[inline]
        fn acos(self) -> Self {
            Real::acos(self)
        }

        #[inline]
        fn atan(self) -> Self {
            Real::atan(self)
        }

        #[inline]
        fn sin_cos(self) -> (Self, Self) {
            Real::sin_cos(self)
        }

//            #[inline]
//            fn exp_m1(self) -> Self {
//                Self::exp_m1(self)
//            }
//
//            #[inline]
//            fn ln_1p(self) -> Self {
//                Self::ln_1p(self)
//            }
//
        #[inline]
        fn sinh(self) -> Self {
            Real::sinh(self)
        }

        #[inline]
        fn cosh(self) -> Self {
            Real::cosh(self)
        }

        #[inline]
        fn tanh(self) -> Self {
            Real::tanh(self)
        }

        #[inline]
        fn asinh(self) -> Self {
            Real::asinh(self)
        }

        #[inline]
        fn acosh(self) -> Self {
            Real::acosh(self)
        }

        #[inline]
        fn atanh(self) -> Self {
            Real::atanh(self)
        }

        #[inline]
        fn is_finite(&self) -> bool {
            crate::Infinite::is_finite(*self)
        }
    }
}

mod simba_impl_finite {
    use simba::scalar::{RealField, SubsetOf};
    use num_traits::{Zero, One};
    use crate::Finite;

    use crate::Real;

    impl<T: simba::simd::PrimitiveSimdValue> simba::simd::PrimitiveSimdValue for crate::Finite<T>
    where
        T: simba::simd::SimdValue
            + crate::Primitive
            + crate::Infinite
            + crate::IntrinsicOrd
            + crate::Nan
            + crate::Real
    {}

    impl<T> simba::simd::SimdValue for crate::Finite<T>
    where
        T: simba::simd::SimdValue
            + crate::Primitive
            + crate::Infinite
            + crate::IntrinsicOrd
            + crate::Nan
            + crate::Real
    {
        type Element = Self;
        type SimdBool = bool;

        #[inline(always)]
        fn lanes() -> usize {
            1
        }

        #[inline(always)]
        fn splat(val: Self::Element) -> Self {
            val
        }

        #[inline(always)]
        fn extract(&self, _: usize) -> Self::Element {
            *self
        }

        #[inline(always)]
        unsafe fn extract_unchecked(&self, _: usize) -> Self::Element {
            *self
        }

        #[inline(always)]
        fn replace(&mut self, _: usize, val: Self::Element) {
            *self = val
        }

        #[inline(always)]
        unsafe fn replace_unchecked(&mut self, _: usize, val: Self::Element) {
            *self = val
        }

        #[inline(always)]
        fn select(self, cond: Self::SimdBool, other: Self) -> Self {
            if cond {
                self
            } else {
                other
            }
        }
    }

    impl<T> simba::scalar::RealField for crate::Finite<T>
    where
        f64: simba::scalar::SubsetOf<T>,
        T: simba::scalar::RealField
            + crate::Primitive
            + crate::Infinite
            + crate::IntrinsicOrd
            + crate::Nan
            + crate::Real
    {
        #[inline]
        fn is_sign_positive(&self) -> bool {
            crate::Encoding::is_sign_positive(self.clone())
        }

        #[inline]
        fn is_sign_negative(&self) -> bool {
            crate::Encoding::is_sign_negative(self.clone())
        }

        #[inline(always)]
        fn copysign(self, to: Self) -> Self {
            // let signbit = (-T::zero()).to_bits();
            // Self::from_bits((signbit & self.to_bits()) | ((!signbit) & to.to_bits()))
            todo!()
        }

        #[inline]
        fn max(self, other: Self) -> Self {
            std::cmp::Ord::max(self, other)
        }

        #[inline]
        fn min(self, other: Self) -> Self {
            std::cmp::Ord::min(self, other)
        }

        #[inline]
        fn clamp(self, min: Self, max: Self) -> Self {
            if self < min {
                min
            } else if self > max {
                max
            } else {
                self
            }
        }

        #[inline]
        fn atan2(self, other: Self) -> Self {
            crate::Real::atan2(self, other)
        }

        /// Archimedes' constant.
        #[inline]
        fn pi() -> Self {
            crate::Real::PI
        }

        /// 2.0 * pi.
        #[inline]
        fn two_pi() -> Self {
            crate::Real::TAU
        }

        /// pi / 2.0.
        #[inline]
        fn frac_pi_2() -> Self {
            crate::Real::FRAC_PI_2
        }

        /// pi / 3.0.
        #[inline]
        fn frac_pi_3() -> Self {
            crate::Real::FRAC_PI_3
        }

        /// pi / 4.0.
        #[inline]
        fn frac_pi_4() -> Self {
            crate::Real::FRAC_PI_4
        }

        /// pi / 6.0.
        #[inline]
        fn frac_pi_6() -> Self {
            crate::Real::FRAC_PI_6
        }

        /// pi / 8.0.
        #[inline]
        fn frac_pi_8() -> Self {
            crate::Real::FRAC_PI_8
        }

        /// 1.0 / pi.
        #[inline]
        fn frac_1_pi() -> Self {
            crate::Real::FRAC_1_PI
        }

        /// 2.0 / pi.
        #[inline]
        fn frac_2_pi() -> Self {
            crate::Real::FRAC_2_PI
        }

        /// 2.0 / sqrt(pi).
        #[inline]
        fn frac_2_sqrt_pi() -> Self {
            crate::Real::FRAC_2_SQRT_PI
        }


        /// Euler's number.
        #[inline]
        fn e() -> Self {
            crate::Real::E
        }

        /// log2(e).
        #[inline]
        fn log2_e() -> Self {
            crate::Real::LOG2_E
        }

        /// log10(e).
        #[inline]
        fn log10_e() -> Self {
            crate::Real::LOG10_E
        }

        /// ln(2.0).
        #[inline]
        fn ln_2() -> Self {
            crate::Real::LN_2
        }

        /// ln(10.0).
        #[inline]
        fn ln_10() -> Self {
            crate::Real::LN_10
        }

        fn min_value() -> Option<Self> {
            T::min_value().and_then(|t| Finite::new(t).ok())
        }

        fn max_value() -> Option<Self> {
            T::max_value().and_then(|t| Finite::new(t).ok())
        }
    }


    impl<T> simba::scalar::SubsetOf<Self> for crate::Finite<T>
    where
        T: simba::scalar::SubsetOf<T>
            + crate::Primitive
            + crate::Infinite
            + crate::IntrinsicOrd
            + crate::Nan
            + crate::Real
    {
        #[inline]
        fn to_superset(&self) -> Self {
            *self
        }

        #[inline]
        fn from_superset_unchecked(element: &Self) -> Self {
            *element
        }

        #[inline]
        fn is_in_subset(_: &Self) -> bool {
            true
        }
    }

    impl<T> simba::scalar::SupersetOf<f64> for crate::Finite<T>
    where
        f64: simba::scalar::SubsetOf<T>,
        T: crate::Primitive
            + crate::Infinite
            + crate::IntrinsicOrd
            + crate::Nan
            + crate::Real
    {
        #[inline]
        fn to_subset(&self) -> Option<f64> {
            f64::from_superset(&self.into_inner())
        }

        #[inline]
        fn is_in_subset(&self) -> bool {
            <f64 as simba::scalar::SubsetOf::<T>>::is_in_subset(&self.into_inner())
        }

        #[inline]
        fn to_subset_unchecked(&self) -> f64 {
            f64::from_superset_unchecked(&self.into_inner())
        }

        #[inline]
        fn from_subset(element: &f64) -> Self {
            crate::Finite::new(element.to_superset()).unwrap()
        }
    }

    impl<T> simba::scalar::SupersetOf<f32> for crate::Finite<T>
    where
        f32: simba::scalar::SubsetOf<T>,
        T: crate::Primitive
            + crate::Infinite
            + crate::IntrinsicOrd
            + crate::Nan
            + crate::Real
    {
        #[inline]
        fn to_subset(&self) -> Option<f32> {
            f32::from_superset(&self.into_inner())
        }

        #[inline]
        fn is_in_subset(&self) -> bool {
            <f32 as simba::scalar::SubsetOf::<T>>::is_in_subset(&self.into_inner())
        }

        #[inline]
        fn to_subset_unchecked(&self) -> f32 {
            f32::from_superset_unchecked(&self.into_inner())
        }

        #[inline]
        fn from_subset(element: &f32) -> Self {
            crate::Finite::new(element.to_superset()).unwrap()
        }
    }

    impl<T> simba::scalar::SubsetOf<T> for crate::Finite<T>
    where T: crate::Primitive
            + crate::Infinite
            + crate::IntrinsicOrd
            + crate::Nan
            + crate::Real
    {
        #[inline]
        fn to_superset(&self) -> T {
            self.into_inner()
        }

        #[inline]
        fn from_superset_unchecked(element: &T) -> Self {
            Self::new_unchecked(*element)
        }

        #[inline]
        fn is_in_subset(element: &T) -> bool {
            !element.is_nan()
        }
    }

    impl<T> simba::scalar::Field for crate::Finite<T>
    where
        T: simba::scalar::Field
            + crate::Primitive
            + crate::Infinite
            + crate::IntrinsicOrd
            + crate::Nan
            + crate::Real
    {}

    impl<T> simba::scalar::ComplexField for crate::Finite<T>
    where
        f64: simba::scalar::SubsetOf<T>,
        T: simba::scalar::RealField
            + simba::scalar::ComplexField
            + crate::Primitive
            + crate::Infinite
            + crate::IntrinsicOrd
            + crate::Nan
            + crate::Real
    {
        type RealField = Self;

        #[inline]
        fn from_real(re: Self::RealField) -> Self {
            re
        }

        #[inline]
        fn real(self) -> Self::RealField {
            self
        }

        #[inline]
        fn imaginary(self) -> Self::RealField {
            Self::zero()
        }

        #[inline]
        fn norm1(self) -> Self::RealField {
            Self::abs(self)
        }

        #[inline]
        fn modulus(self) -> Self::RealField {
            Self::abs(self)
        }

        #[inline]
        fn modulus_squared(self) -> Self::RealField {
            self * self
        }

        #[inline]
        fn argument(self) -> Self::RealField {
            if self >= Self::zero() {
                Self::zero()
            } else {
                Self::pi()
            }
        }

        #[inline]
        fn to_exp(self) -> (Self, Self) {
            if self >= Self::zero() {
                (self, Self::one())
            } else {
                (-self, -Self::one())
            }
        }

        #[inline]
        fn recip(self) -> Self {
            crate::Real::recip(self)
        }

        #[inline]
        fn conjugate(self) -> Self {
            self
        }

        #[inline]
        fn scale(self, factor: Self::RealField) -> Self {
            self * factor
        }

        #[inline]
        fn unscale(self, factor: Self::RealField) -> Self {
            self / factor
        }

        #[inline]
        fn floor(self) -> Self {
            crate::Real::floor(self)
        }

        #[inline]
        fn ceil(self) -> Self {
            crate::Real::ceil(self)
        }

        #[inline]
        fn round(self) -> Self {
            crate::Real::round(self)
        }

        #[inline]
        fn trunc(self) -> Self {
            crate::Real::trunc(self)
        }

        #[inline]
        fn fract(self) -> Self {
            crate::Real::fract(self)
        }

        #[inline]
        fn abs(self) -> Self {
            crate::Signed::abs(&self)
        }

        #[inline]
        fn signum(self) -> Self {
            crate::Signed::signum(&self)
        }

        #[inline]
        fn mul_add(self, a: Self, b: Self) -> Self {
            crate::Real::mul_add(self, a, b)
        }

        #[cfg(feature = "std")]
        #[inline]
        fn powi(self, n: i32) -> Self {
            crate::Real::powi(self, n)
        }

        #[cfg(not(feature = "std"))]
        #[inline]
        fn powi(self, n: i32) -> Self {
            // FIXME: is there a more accurate solution?
            Self::powf(self, n as Self)
        }

        #[inline]
        fn powf(self, n: Self) -> Self {
            Real::powf(self, n)
        }

        #[inline]
        fn powc(self, n: Self) -> Self {
            // Same as powf.
            Real::powf(self, n)
        }

        #[inline]
        fn sqrt(self) -> Self {
            Real::sqrt(self)
        }

        #[inline]
        fn try_sqrt(self) -> Option<Self> {
            if self >= Self::zero() {
                Some(Real::sqrt(self))
            } else {
                None
            }
        }

        #[inline]
        fn exp(self) -> Self {
            Real::exp(self)
        }

        #[inline]
        fn exp2(self) -> Self {
            Real::exp2(self)
        }


        #[inline]
        fn exp_m1(self) -> Self {
            Real::exp_m1(self)
        }

        #[inline]
        fn ln_1p(self) -> Self {
            Real::ln_1p(self)
        }

        #[inline]
        fn ln(self) -> Self {
            Real::ln(self)
        }

        #[inline]
        fn log(self, base: Self) -> Self {
            Real::log(self, base)
        }

        #[inline]
        fn log2(self) -> Self {
            Real::log2(self)
        }

        #[inline]
        fn log10(self) -> Self {
            Real::log10(self)
        }

        #[inline]
        fn cbrt(self) -> Self {
            Real::cbrt(self)
        }

        #[inline]
        fn hypot(self, other: Self) -> Self::RealField {
            Real::hypot(self, other)
        }

        #[inline]
        fn sin(self) -> Self {
            Real::sin(self)
        }

        #[inline]
        fn cos(self) -> Self {
            Real::cos(self)
        }

        #[inline]
        fn tan(self) -> Self {
            Real::tan(self)
        }

        #[inline]
        fn asin(self) -> Self {
            Real::asin(self)
        }

        #[inline]
        fn acos(self) -> Self {
            Real::acos(self)
        }

        #[inline]
        fn atan(self) -> Self {
            Real::atan(self)
        }

        #[inline]
        fn sin_cos(self) -> (Self, Self) {
            Real::sin_cos(self)
        }

//            #[inline]
//            fn exp_m1(self) -> Self {
//                Self::exp_m1(self)
//            }
//
//            #[inline]
//            fn ln_1p(self) -> Self {
//                Self::ln_1p(self)
//            }
//
        #[inline]
        fn sinh(self) -> Self {
            Real::sinh(self)
        }

        #[inline]
        fn cosh(self) -> Self {
            Real::cosh(self)
        }

        #[inline]
        fn tanh(self) -> Self {
            Real::tanh(self)
        }

        #[inline]
        fn asinh(self) -> Self {
            Real::asinh(self)
        }

        #[inline]
        fn acosh(self) -> Self {
            Real::acosh(self)
        }

        #[inline]
        fn atanh(self) -> Self {
            Real::atanh(self)
        }

        #[inline]
        fn is_finite(&self) -> bool {
            true
        }
    }
}