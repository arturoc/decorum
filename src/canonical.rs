use core::mem;

use crate::{Encoding, Nan};

const SIGN_MASK: u64 = 0x8000_0000_0000_0000;
const EXPONENT_MASK: u64 = 0x7ff0_0000_0000_0000;
const MANTISSA_MASK: u64 = 0x000f_ffff_ffff_ffff;

const CANONICAL_NAN_BITS: u64 = 0x7ff8_0000_0000_0000;
const CANONICAL_ZERO_BITS: u64 = 0x0;

pub trait ToCanonicalBits: Copy + Sized {
    fn to_canonical_bits(self) -> u64;
}

impl<T> ToCanonicalBits for T
where
    T: Encoding + Nan,
{
    fn to_canonical_bits(self) -> u64 {
        if self.is_nan() {
            CANONICAL_NAN_BITS
        }
        else {
            let (mantissa, exponent, sign) = self.integer_decode();
            if mantissa == 0 {
                CANONICAL_ZERO_BITS
            }
            else {
                let exponent = u64::from(unsafe { mem::transmute::<i16, u16>(exponent) });
                let sign = if sign > 0 { 1u64 } else { 0u64 };
                (mantissa & MANTISSA_MASK)
                    | ((exponent << 52) & EXPONENT_MASK)
                    | ((sign << 63) & SIGN_MASK)
            }
        }
    }
}
