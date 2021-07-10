#![no_std]

extern crate alloc;
use alloc::vec::Vec;

use core::ops::{Bound, RangeBounds};

mod sieve;

pub struct PrimeIterator {
    start: u64,
    end: u64,
}

impl Iterator for PrimeIterator {
    type Item = u64;
    fn next(&mut self) -> Option<Self::Item> {
        None
    }
}

pub fn get_primes_up_to(n: u64) -> Vec<u64> {
    get_primes_in_range(..=n)
}

pub fn iterate_primes_up_to(n: u64) -> PrimeIterator {
    iterate_primes_in_range(..=n)
}

pub fn get_primes_in_range<R>(range: R) -> Vec<u64>
where
    R: RangeBounds<u64>,
{
    let start = match range.start_bound() {
        Bound::Included(s) => *s,
        Bound::Excluded(s) => s.saturating_add(1),
        Bound::Unbounded => 0,
    };
    let end = match range.end_bound() {
        Bound::Included(s) => s.saturating_add(1),
        Bound::Excluded(s) => *s,
        Bound::Unbounded => u64::MAX,
    };
    todo!()
}

pub fn iterate_primes_in_range<R>(range: R) -> PrimeIterator
where
    R: RangeBounds<u64>,
{
    let start = match range.start_bound() {
        Bound::Included(s) => *s,
        Bound::Excluded(s) => s.saturating_add(1),
        Bound::Unbounded => 0,
    };
    let end = match range.end_bound() {
        Bound::Included(s) => s.saturating_add(1),
        Bound::Excluded(s) => *s,
        Bound::Unbounded => u64::MAX,
    };
    todo!()
}
