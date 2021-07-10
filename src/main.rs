#![warn(unsafe_op_in_unsafe_fn)]

use std::{io::BufRead, time::Instant};

fn main() {
    std::io::stdin()
        .lock()
        .read_line(&mut "".to_string())
        .unwrap();

    for e in 3..=10 {
        let n = 10_usize.pow(e);
        println!();
        println!("n = {}", n);

        let start = Instant::now();
        let correct = sqrt_sieve(n);
        println!("Sqrt sieve: {} ms", start.elapsed().as_millis());
        // let start = Instant::now();
        // let correct = segmented_wheel30_sieve(n);
        // println!("Segmented wheel30 sieve: {} ms", start.elapsed().as_millis());
        println!(
            "Segmented wheel30 sieve: {} ms",
            time_it(n, segmented_wheel30_sieve, &correct)
        );
    }
}

fn time_it(n: usize, f: fn(usize) -> Vec<u64>, correct: &[u64]) -> u128 {
    let start = Instant::now();
    let res = f(n);
    let time = start.elapsed().as_millis();
    for i in 0..correct.len() {
        if correct[i] != res[i] {
            dbg!(i, correct[i], res[i]);
            panic!()
        }
    }
    assert_eq!(&res, correct);
    time
}

fn sqrt_sieve(n: usize) -> Vec<u64> {
    let mut is_prime = vec![true; n + 1];
    // is_prime[0] = false;
    // is_prime[1] = false;

    (2..=n)
        .filter_map(|p| {
            is_prime[p].then(|| {
                for multiple in (p * p..=n).step_by(p) {
                    is_prime[multiple] = false;
                }
                p as u64
            })
        })
        .collect()

    // is_prime
    //     .into_iter()
    //     .enumerate()
    //     .filter_map(|(p, b)| if b { Some(p as u64) } else { None })
    //     .collect()
}

const SMALL_PRIMES: [u64; 410] = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
    101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193,
    197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307,
    311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421,
    431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547,
    557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659,
    661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797,
    809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929,
    937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033, 1039,
    1049, 1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153,
    1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249, 1259, 1277, 1279,
    1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 1381, 1399, 1409,
    1423, 1427, 1429, 1433, 1439, 1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499,
    1511, 1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 1597, 1601, 1607, 1609, 1613,
    1619, 1621, 1627, 1637, 1657, 1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, 1741,
    1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, 1823, 1831, 1847, 1861, 1867, 1871, 1873,
    1877, 1879, 1889, 1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987, 1993, 1997, 1999,
    2003, 2011, 2017, 2027, 2029, 2039, 2053, 2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113,
    2129, 2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213, 2221, 2237, 2239, 2243, 2251,
    2267, 2269, 2273, 2281, 2287, 2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, 2371,
    2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, 2437, 2441, 2447, 2459, 2467, 2473, 2477,
    2503, 2521, 2531, 2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, 2621, 2633, 2647,
    2657, 2659, 2663, 2671, 2677, 2683, 2687, 2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731,
    2741, 2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819,
];
const SMALL_SIEVING_PRIMES: [&[u64]; 8] = [
    &SMALL_SIEVING_PRIMES_1,
    &SMALL_SIEVING_PRIMES_7,
    &SMALL_SIEVING_PRIMES_11,
    &SMALL_SIEVING_PRIMES_13,
    &SMALL_SIEVING_PRIMES_17,
    &SMALL_SIEVING_PRIMES_19,
    &SMALL_SIEVING_PRIMES_23,
    &SMALL_SIEVING_PRIMES_29,
];
const SMALL_SIEVING_PRIMES_1: [u64; 48] = [
    1, 2, 5, 6, 7, 8, 9, 11, 14, 18, 19, 20, 21, 22, 23, 25, 27, 33, 34, 35, 39, 40, 41, 43, 44,
    46, 49, 51, 54, 58, 60, 61, 62, 65, 67, 71, 72, 74, 75, 76, 77, 78, 79, 84, 85, 89, 91, 93,
];
const SMALL_SIEVING_PRIMES_7: [u64; 52] = [
    0, 1, 2, 3, 4, 5, 9, 10, 11, 12, 13, 15, 16, 18, 19, 20, 24, 25, 26, 29, 30, 31, 32, 33, 36,
    37, 41, 43, 44, 48, 52, 53, 54, 55, 58, 59, 62, 66, 67, 71, 76, 78, 79, 81, 82, 85, 87, 88, 89,
    90, 92, 93,
];
const SMALL_SIEVING_PRIMES_11: [u64; 52] = [
    0, 1, 2, 3, 4, 6, 8, 9, 10, 13, 14, 15, 16, 17, 21, 23, 25, 27, 29, 30, 31, 32, 34, 35, 36, 38,
    39, 43, 45, 48, 49, 50, 52, 53, 57, 60, 62, 63, 64, 69, 70, 71, 78, 79, 80, 81, 84, 86, 87, 90,
    91, 93,
];
const SMALL_SIEVING_PRIMES_13: [u64; 52] = [
    0, 1, 2, 3, 5, 6, 7, 9, 10, 12, 14, 15, 17, 20, 21, 22, 24, 27, 28, 29, 34, 35, 36, 37, 38, 40,
    43, 47, 48, 49, 51, 55, 56, 57, 58, 59, 62, 64, 66, 68, 69, 70, 71, 73, 76, 79, 82, 83, 86, 89,
    90, 93,
];
const SMALL_SIEVING_PRIMES_17: [u64; 52] = [
    0, 1, 3, 4, 5, 6, 7, 8, 10, 11, 15, 18, 19, 20, 21, 22, 26, 27, 28, 29, 31, 32, 36, 39, 40, 42,
    43, 45, 47, 49, 53, 54, 55, 56, 59, 61, 62, 63, 66, 67, 69, 73, 74, 75, 76, 78, 80, 81, 82, 88,
    89, 92,
];
const SMALL_SIEVING_PRIMES_19: [u64; 48] = [
    0, 2, 3, 4, 6, 7, 11, 12, 13, 14, 16, 20, 23, 24, 25, 27, 28, 30, 33, 34, 35, 37, 41, 42, 46,
    47, 48, 49, 51, 52, 53, 55, 56, 58, 59, 62, 66, 67, 69, 72, 74, 75, 79, 84, 88, 89, 90, 91,
];
const SMALL_SIEVING_PRIMES_23: [u64; 52] = [
    0, 1, 2, 3, 5, 7, 8, 9, 11, 12, 14, 16, 18, 19, 21, 22, 24, 25, 28, 31, 32, 33, 36, 38, 39, 40,
    42, 45, 47, 49, 50, 51, 52, 53, 57, 60, 63, 65, 66, 68, 71, 73, 74, 75, 77, 79, 80, 84, 87, 88,
    89, 91,
];
const SMALL_SIEVING_PRIMES_29: [u64; 51] = [
    0, 1, 2, 4, 5, 7, 8, 11, 12, 13, 14, 15, 16, 18, 19, 21, 23, 26, 27, 30, 33, 34, 36, 40, 41,
    42, 43, 46, 47, 49, 51, 53, 56, 62, 64, 65, 67, 68, 69, 70, 76, 77, 79, 81, 84, 85, 86, 89, 90,
    92, 93,
];
const MAX_PRECOMPUTED: usize = 2833;
const SKIP: usize = MAX_PRECOMPUTED / 30;
const SEGMENT_SIZE: usize = 2_usize.pow(18); // should be a big enough power of two
const REMS: [usize; 8] = [1, 7, 11, 13, 17, 19, 23, 29];
const REM_MUL_DIV: [[(usize, u8); 8]; 8] = {
    let mut arr = [[(0, 0); 8]; 8];
    let mut i = 0;
    while i < 8 {
        let mut j = 0;
        while j < 8 {
            let prod = REMS[i] * REMS[j];
            arr[i][j] = (prod / 30, 1 << MUL_REM[i][j]);
            j += 1;
        }
        i += 1;
    }
    arr
};
const MUL_REM: [[usize; 8]; 8] = {
    let mut arr = [[0; 8]; 8];
    let mut i = 0;
    while i < 8 {
        let mut j = 0;
        while j < 8 {
            arr[i][j] = rem_to_index((REMS[i] * REMS[j]) as u64 % 30);
            j += 1;
        }
        i += 1;
    }
    arr
};

const fn rem_to_index(n: u64) -> usize {
    match n {
        1 => 0,
        7 => 1,
        11 => 2,
        13 => 3,
        17 => 4,
        19 => 5,
        23 => 6,
        29 => 7,
        _ => 9,
    }
}

#[inline(never)]
fn segmented_wheel30_sieve(n: usize) -> Vec<u64> {
    if n < MAX_PRECOMPUTED {
        let prime_count = count_less_than(n, &SMALL_PRIMES);
        return SMALL_PRIMES[..prime_count].to_vec();
    }

    // TODO: add sieving primes if needed
    let (mut primes, mut is_prime_buffer, mut sieving_primes) = sieve_first_segment(n);

    for cur_segment in 1..=n / (SEGMENT_SIZE * 30) {
        is_prime_buffer.fill(u8::MAX);
        let start = SEGMENT_SIZE * cur_segment;
        sieve_segment(
            start as u64,
            &mut is_prime_buffer,
            &[
                &sieving_primes[0],
                &sieving_primes[1],
                &sieving_primes[2],
                &sieving_primes[3],
                &sieving_primes[4],
                &sieving_primes[5],
                &sieving_primes[6],
                &sieving_primes[7],
            ],
        );
        unsafe {
            append_primes(&mut primes, &is_prime_buffer, start);
        }
    }

    let prime_count = count_less_than(n, &primes);
    primes.truncate(prime_count);

    primes
}

fn count_less_than(n: usize, primes: &[u64]) -> usize {
    match primes.binary_search(&(n as u64)) {
        Ok(i) => i + 1,
        Err(i) => i,
    }
}

fn prime_count_upper_bound(n: usize) -> usize {
    ((n + 30) as f64 / (((n + 30) as f64).ln() - 2.0)) as usize
}

fn sieve_first_segment(n: usize) -> (Vec<u64>, Vec<u8>, [Vec<u64>; 8]) {
    let len = SEGMENT_SIZE.min(n / 30 + 1);
    let mut is_prime = vec![u8::MAX; len];
    let extra = if len == SEGMENT_SIZE {
        8 * SEGMENT_SIZE
    } else {
        0
    };
    let mut primes = Vec::with_capacity(prime_count_upper_bound(n) + extra);
    primes.extend_from_slice(&SMALL_PRIMES);

    let mut sieving_primes = [
        SMALL_SIEVING_PRIMES_1.to_vec(),
        SMALL_SIEVING_PRIMES_7.to_vec(),
        SMALL_SIEVING_PRIMES_11.to_vec(),
        SMALL_SIEVING_PRIMES_13.to_vec(),
        SMALL_SIEVING_PRIMES_17.to_vec(),
        SMALL_SIEVING_PRIMES_19.to_vec(),
        SMALL_SIEVING_PRIMES_23.to_vec(),
        SMALL_SIEVING_PRIMES_29.to_vec(),
    ];

    sieve_segment(
        SKIP as u64,
        &mut is_prime[..len - SKIP],
        &SMALL_SIEVING_PRIMES,
    );

    unsafe {
        append_primes(&mut primes, &is_prime[..len - SKIP], SKIP);
        if MAX_PRECOMPUTED * MAX_PRECOMPUTED < n {
            append_sieving_primes(&mut sieving_primes, &is_prime[..len - SKIP], SKIP);
        }
    }

    (primes, is_prime, sieving_primes)
}

fn sieve_segment(start: u64, is_prime: &mut [u8], sieving_primes: &[&[u64]; 8]) {
    sieve_segment_with_rem::<0>(start, is_prime, sieving_primes[0]);
    sieve_segment_with_rem::<1>(start, is_prime, sieving_primes[1]);
    sieve_segment_with_rem::<2>(start, is_prime, sieving_primes[2]);
    sieve_segment_with_rem::<3>(start, is_prime, sieving_primes[3]);
    sieve_segment_with_rem::<4>(start, is_prime, sieving_primes[4]);
    sieve_segment_with_rem::<5>(start, is_prime, sieving_primes[5]);
    sieve_segment_with_rem::<6>(start, is_prime, sieving_primes[6]);
    sieve_segment_with_rem::<7>(start, is_prime, sieving_primes[7]);
}

#[inline(never)]
fn sieve_segment_with_rem<const R_INDEX: usize>(
    start: u64,
    is_prime: &mut [u8],
    sieving_primes: &[u64],
) {
    let len = is_prime.len() as u64;
    'prime_loop: for p in sieving_primes
        .iter()
        .copied()
        .take_while(|p| 30 * p * p <= start + len)
    {
        let q = p as usize;
        let p = (30 * q + REMS[R_INDEX]) as u64;
        let seg_rem = (start % p) as usize;
        let offset_mask_iter = (0..8).map(|i| {
            (
                q * REMS[i] + REM_MUL_DIV[R_INDEX][i].0,
                !REM_MUL_DIV[R_INDEX][i].1,
            )
        });

        for (offset, mask) in offset_mask_iter
            .clone()
            .skip_while(|(offset, _)| *offset < seg_rem)
        {
            if let Some(b) = is_prime.get_mut(offset - seg_rem) {
                *b &= mask;
            } else {
                continue 'prime_loop;
            }
        }

        // let mut i_base = p as usize - seg_rem;
        // let count = ((len as usize).saturating_sub(i_base + REM_MUL_DIV[R_INDEX][7].0 + 29 * q)) / p as usize;
        // if count != 0 {
        //     let mut is_prime_ptr = unsafe { is_prime.as_mut_ptr().add(i_base) };
        //     for _ in 0..count {
        //         for (offset, mask) in offset_mask_iter.clone() {
        //             unsafe {
        //                 debug_assert!(is_prime_ptr < is_prime.as_mut_ptr().add(is_prime.len()));
        //                 *is_prime_ptr.add(offset) &= mask;
        //             }
        //         }
        //         is_prime_ptr = unsafe { is_prime_ptr.add(p as usize) };
        //     }
        //
        //     i_base += count * p as usize;
        // }

        let mut i_base = p as usize - seg_rem;
        let mut is_prime_ptr = unsafe { is_prime.as_mut_ptr().add(i_base) };
        let max_ptr = unsafe {
            is_prime
                .as_mut_ptr()
                .add(is_prime.len())
                .sub(offset_mask_iter.clone().last().unwrap().0)
        };
        while is_prime_ptr < max_ptr {
            for (offset, mask) in offset_mask_iter.clone() {
                unsafe {
                    debug_assert!(is_prime_ptr < is_prime.as_mut_ptr().add(is_prime.len()));
                    *is_prime_ptr.add(offset) &= mask;
                }
            }
            is_prime_ptr = unsafe { is_prime_ptr.add(p as usize) };
        }
        i_base = is_prime_ptr as usize - is_prime.as_mut_ptr() as usize;

        // let mut i_base = p as usize - seg_rem;
        // while i_base + REM_MUL_DIV[R_INDEX][7].0 + 29 * q < is_prime.len() {
        //     for (offset, mask) in offset_mask_iter.clone() {
        //         unsafe {
        //             *is_prime.get_unchecked_mut(i_base + offset) &= mask;
        //         }
        //     }
        //     i_base += p as usize;
        // }

        // REMS[i] <= 29
        // offset
        // = q * REMS[i] + REM_MUL_DIV[R_INDEX][i].0
        // = floor(p / 30) * REMS[i] + floor( ( (p % 30) * REMS[i] ) / 30 )
        // <= floor(p / 30) * REMS[i] + ((p % 30) * REMS[i]) / 30
        // = REMS[i] * (floor(p / 30) + (p % 30) / 30)
        // = (REMS[i] / 30) * (30 * floor(p / 30) + p % 30)
        // = (REMS[i] / 30) * p
        // < p

        // let mut i_base = p as usize - seg_rem;
        // for chunk in is_prime[i_base..].chunks_exact_mut(p as usize) {
        //     for (offset, mask) in offset_mask_iter.clone() {
        //         unsafe {
        //             *chunk.get_unchecked_mut(offset) &= mask;
        //         }
        //     }
        // }

        // i_base += is_prime[i_base..].chunks_exact_mut(p as usize).count() * p as usize;

        let mut offset_mask = offset_mask_iter.clone().next().unwrap();
        let mut j = 0;
        while let Some(b) = is_prime.get_mut(i_base + offset_mask.0) {
            *b &= offset_mask.1;
            j += 1;
            i_base += (j / 8) * (p as usize);
            j %= 8;
            offset_mask = offset_mask_iter.clone().nth(j).unwrap();
        }
    }
}

unsafe fn append_primes(primes: &mut Vec<u64>, is_prime: &[u8], start: usize) {
    let mut cur_num = 30 * start;
    let mut end = unsafe { primes.as_mut_ptr().add(primes.len()) };
    for byte in is_prime.iter().copied() {
        for (shift, rem) in REMS.iter().copied().enumerate() {
            debug_assert!((end as usize - primes.as_ptr() as usize) / 8 < primes.capacity());
            unsafe {
                *end = (cur_num + rem) as u64;
                end = end.add(1 & (byte >> shift) as usize);
            }
        }
        cur_num += 30;
    }

    assert!((end as usize - primes.as_ptr() as usize) / 8 < primes.capacity());
    unsafe {
        primes.set_len((end as usize - primes.as_ptr() as usize) / 8);
    }
}

unsafe fn append_sieving_primes(sieving_primes: &mut [Vec<u64>; 8], is_prime: &[u8], start: usize) {
    let mut cur_num = start as u64;
    for byte in is_prime.iter().copied() {
        for shift in 0..8 {
            if byte & (1 << shift) != 0 {
                sieving_primes[shift].push(cur_num);
            }
        }
        cur_num += 1;
    }
}
