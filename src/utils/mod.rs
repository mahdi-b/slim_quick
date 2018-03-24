#![allow(non_snake_case)]
#![allow(unused_variables)]
#![allow(dead_code)]
extern crate ndarray;

use rand::{thread_rng, SeedableRng, sample, StdRng};

use self::ndarray::prelude::*;

//////////////////////////////////////////////////////
//  MARK DNA
//////////////////////////////////////////////////////
pub enum DNA {A,C,G,T} //0,1,2,3

static PRIMES_257: [u16; 257] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 1597, 1601, 1607, 1609, 1613, 1619, 1621];
static LOG_PRIMES_257: [u16; 257] = [0; 257];

//////////////////////////////////////////////////////
//  MARK Sequence Utils
//////////////////////////////////////////////////////


pub fn convert_byte_DNA(byte_val: &u8) -> Option<usize> {
    match *byte_val {
        65 => Some(DNA::A as usize),
        67 => Some(DNA::C as usize),
        71 => Some(DNA::G as usize),
        84 => Some(DNA::T as usize),

        // TODO: handle cases where N is randomly selected as either an base
        // for now, consider N as DNA::A
        78 => Some(DNA::A as usize),
        _ => None
    }
}


// pre-computes the random samples that will be used for subset.
// each subset has a samnple of size subset_size.
pub fn get_random_samples(kmer_size:usize, subset_size: usize, nb_subsets: usize, seed:Option<&[usize]>) -> Vec<Vec<u16>> {

    let mut bands = vec![vec![0u16; subset_size ]; nb_subsets];

    if let Some(s) = seed {
        let mut rng: StdRng = SeedableRng::from_seed(s);
        for i in 0..nb_subsets{
            bands[i] = sample(&mut rng, 0..(4u16.pow(kmer_size as u32)), subset_size);
        }
    }else{
        let mut rng = thread_rng();
        for i in 0..nb_subsets{
            bands[i] = sample(&mut rng, 0..(4u16.pow(kmer_size as u32)), subset_size);
        }
    }
    return bands;
}

pub fn build_filters_matrix(nb_kmers: usize, nb_subsets:usize, bands_vec:Vec<f64>){

    let mut log_p: Vec<f64> = Vec::new();

    for (i, val) in  PRIMES_257.iter().enumerate()  {
        if i == nb_kmers{ break;}

        log_p.push((*val as f64).log10());

    }

    let filters_mtx = Array::from_shape_vec((nb_kmers,nb_subsets), bands_vec).unwrap();


    let logPrimesVec = 	Array::from_shape_vec((nb_kmers, 1), log_p[..nb_kmers].to_vec()).unwrap();

    let log_p_mtx = filters_mtx * logPrimesVec;


}


pub struct Utils {
    pub log_primes: Vec<f64>,
}

//
impl Utils {
    pub fn new() -> Utils {
        let mut log_p: Vec<f64> = Vec::new();
        for (i, val) in  PRIMES_257.iter().enumerate()  {
            log_p.push((*val as f64).log10());
        }
        Utils {log_primes: log_p}
    }
}


////////////////////////////////////////////////////////
//// Mark Tests
////////////////////////////////////////////////////////

#[cfg(test)]
mod test {

    mod test_get_random_samples{
        use super::super::*;
        #[test]
        fn test_get_random_sample() {
            let (kmer_size, subset_size, nb_subsets) = (2, 3,6);
            let random_samples = get_random_samples(kmer_size, subset_size, nb_subsets, Some(&[1,2,3,4]));

            println!("{:?}", random_samples);
            assert_eq!(random_samples, [[15, 10, 3], [0, 1, 12], [3, 6, 5], [12, 1, 9], [5, 15, 12], [0, 13, 9]])
        }
    }


    //////////////////////////////////////////////////////
    // MARK convert_byte_DNA tests
    //////////////////////////////////////////////////////
    mod test_convert_by_DNA {
        use super::super::*;
        #[test]
        fn returns_A() {
            assert_eq!(Some(0), convert_byte_DNA(&65))
        }

        #[test]
        fn returns_C() {
            assert_eq!(Some(1), convert_byte_DNA(&67))
        }

        #[test]
        fn returns_G() {
            assert_eq!(Some(2), convert_byte_DNA(&71))
        }

        #[test]
        fn returns_T() {
            assert_eq!(Some(3), convert_byte_DNA(&84))
        }

        #[test]
        fn returns_other() {
            assert_eq!(None, convert_byte_DNA(&0))
        }
    }
    //////////////////////////////////////////////////////
    // MARK utils tests
    //////////////////////////////////////////////////////

    mod utils_test{
        use super::super::*;

        #[test]
        fn test_new(){
            let utils = Utils::new();
            // TODO add appropriate test here for the Utils: strucut
            assert_eq!(1,1)
        }

    }

}