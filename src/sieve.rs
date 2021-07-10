extern crate alloc;
use alloc::vec::Vec;

struct PrimeSieve {
    small_segment_size: usize,
    big_segment_size: usize,
    small_prime_threshold: u64,
    medium_prime_threshold: u64,
}

impl PrimeSieve {
    fn sieve_range(&self, start: u64, end: u64) -> SievedPrimes {
        if start == end {
            return todo!();
        }
        assert!(start < end);
        
        let (small_sieving_primes, medium_sieving_primes) = self.get_sieving_primes(end);
        
        
        
        todo!();
    }
    
    fn get_sieving_primes(&self, end: u64) -> (Vec<u64>, [Vec<u64>; 8]) {
        todo!()
    }
    
    fn sieve_range_with_buffer(&self, start: u64, end: u64, buffer: &mut [u8]) {
        if start < self.big_segment_size as u64 {
            
        }
    }   

    fn sieve_segment(&self, start: u64, is_prime: &mut [u8], sieving_primes: &[Vec<u64>; 8]) {

    }
}

impl Default for PrimeSieve {
    fn default() -> PrimeSieve {
        PrimeSieve {
            small_segment_size: 2_usize.pow(16),
            big_segment_size: 2_usize.pow(18),
            small_prime_threshold: 2048,
            medium_prime_threshold: 100_000,
        }
    }
}

struct SievedPrimes {
    start: u64, 
    end: u64,
    buffer: Vec<u8>,
}