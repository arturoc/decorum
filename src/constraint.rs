//! Constraints on the members of floating-point values that proxy types may
//! represent.

use core::convert::Infallible;
use core::fmt::Debug;
use core::marker::PhantomData;
#[cfg(feature = "std")]
use thiserror::Error;

use crate::{Float, Primitive};

#[cfg_attr(feature = "std", derive(Error))]
#[cfg_attr(feature = "std", error("floating-point constraint violated"))]
#[derive(Clone, Copy, Debug)]
pub struct ConstraintViolation;

pub trait ExpectConstrained<T>: Sized {
    fn expect_constrained(self) -> T;
}

impl<T, E> ExpectConstrained<T> for Result<T, E>
where
    E: Debug,
{
    #[cfg(not(feature = "std"))]
    fn expect_constrained(self) -> T {
        self.expect("floating-point constraint violated")
    }

    #[cfg(feature = "std")]
    fn expect_constrained(self) -> T {
        // In `std` environments, `ConstraintViolation` implements `Error` and
        // an appropriate error message is displayed.
        self.unwrap()
    }
}

pub enum RealSet {}
pub enum InfiniteSet {}
pub enum NanSet {}

pub trait Member<T> {}

pub trait SupersetOf<P> {}

pub trait SubsetOf<P> {}

impl<P, Q> SubsetOf<Q> for P where Q: SupersetOf<P> {}

/// Describes constraints on the set of floating-point values that a proxy type
/// may represent.
///
/// This trait expresses a constraint by filter-mapping values. Note that
/// constraints require `Member<RealSet>`, meaning that the set of real numbers
/// must always be supported and is implied.
pub trait Constraint<T>: Member<RealSet>
where
    T: Float + Primitive,
{
    type Error: Debug;

    /// Filter-maps a primitive floating-point value based on some constraints.
    ///
    /// Returns `None` for values that cannot satify constraints.
    fn filter_map(inner: T) -> Result<T, Self::Error>;
}

#[derive(Debug)]
pub struct UnitConstraint<T>
where
    T: Float + Primitive,
{
    phantom: PhantomData<fn() -> T>,
}

// TODO: Should implementations map values like zero and `NaN` to canonical
//       forms?
impl<T> Constraint<T> for UnitConstraint<T>
where
    T: Float + Primitive,
{
    type Error = Infallible;

    fn filter_map(inner: T) -> Result<T, Self::Error> {
        Ok(inner)
    }
}

impl<T> Member<InfiniteSet> for UnitConstraint<T> where T: Float + Primitive {}

impl<T> Member<NanSet> for UnitConstraint<T> where T: Float + Primitive {}

impl<T> Member<RealSet> for UnitConstraint<T> where T: Float + Primitive {}

impl<T> SupersetOf<FiniteConstraint<T>> for UnitConstraint<T> where T: Float + Primitive {}

impl<T> SupersetOf<NotNanConstraint<T>> for UnitConstraint<T> where T: Float + Primitive {}

/// Disallows `NaN`s.
#[derive(Debug)]
pub struct NotNanConstraint<T>
where
    T: Float + Primitive,
{
    phantom: PhantomData<fn() -> T>,
}

impl<T> Constraint<T> for NotNanConstraint<T>
where
    T: Float + Primitive,
{
    type Error = ConstraintViolation;

    fn filter_map(inner: T) -> Result<T, Self::Error> {
        if inner.is_nan() {
            Err(ConstraintViolation)
        }
        else {
            Ok(inner)
        }
    }
}

impl<T> Member<InfiniteSet> for NotNanConstraint<T> where T: Float + Primitive {}

impl<T> Member<RealSet> for NotNanConstraint<T> where T: Float + Primitive {}

impl<T> SupersetOf<FiniteConstraint<T>> for NotNanConstraint<T> where T: Float + Primitive {}

/// Disallows `NaN`s and infinities.
#[derive(Debug)]
pub struct FiniteConstraint<T>
where
    T: Float + Primitive,
{
    phantom: PhantomData<fn() -> T>,
}

impl<T> Constraint<T> for FiniteConstraint<T>
where
    T: Float + Primitive,
{
    type Error = ConstraintViolation;

    fn filter_map(inner: T) -> Result<T, Self::Error> {
        if inner.is_nan() || inner.is_infinite() {
            Err(ConstraintViolation)
        }
        else {
            Ok(inner)
        }
    }
}

impl<T> Member<RealSet> for FiniteConstraint<T> where T: Float + Primitive {}
