#![allow(unused_variables)]
#![allow(unused_mut)]
#![allow(unused_imports)]
//TODO: replace the word band with subset throughout the code
use std::collections::btree_map::BTreeMap;

extern crate slim_quick;

#[macro_use]
extern crate ndarray;
use ndarray::prelude::*;


#[macro_use]
extern crate lazy_static;

extern crate rand;

//mod data_structures;
//mod utils;
use slim_quick::data_structures;
use slim_quick::utils;



use std::io;
use std::collections::HashMap;
use std::collections::HashSet;


extern crate bio;
use bio::io::fasta;
use bio::io::fasta::Record;
use bio::data_structures::bitenc::BitEnc;


extern crate multiset;
use multiset::HashMultiSet;

extern crate time;
use time::PreciseTime;


extern crate bloomfilter;
use bloomfilter::Bloom;

use rand::{thread_rng, SeedableRng, sample, StdRng};

//
//fn getFilter(bands: Vec<Vec<u16>>, k=u8){
//
//    let nb_kemrs = 4u8.pow(k);
//
//    let  nb_subsets = bands.len();
//    filter_mtx = Array::zeros((nb_kemrs, nb_subsets));
//
//    //TO DO: find a better way to do this?
//    for (col, rows) in  bands.iter().enumerate:
//        filter_mtx[[rows, col]] = 1;
//    return np.row_stack((LOG_PRIMES[:nb_subsets].reshape(1,-1)[0], filter_mtx))
//
//
//}



fn logPrimesVect(nb_kmers: usize, ) -> Vec<f64>{

    static PRIMES_257: [u16; 257] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 1597, 1601, 1607, 1609, 1613, 1619, 1621];


    let mut log_p: Vec<f64> = Vec::new();
    for (i, val) in  PRIMES_257.iter().enumerate()  {
        if i == nb_kmers{ break;}

        log_p.push((*val as f64).log10());

    }
    log_p

}

pub fn get_random_samples(nb_kmers: usize, subset_size: usize, nb_subsets: usize, seed:Option<&[usize]>) -> Vec<f64> {

    let mut bands_vec = vec![0.0f64; (nb_kmers * nb_subsets) as usize];
    let mut rng = thread_rng();
    for i in 0..nb_subsets {
        for val in sample(&mut rng, 0..nb_kmers, subset_size){
            //println!("subseet:{} val:{}", i, val);
            let index = nb_subsets * val + i;
            bands_vec[ index as usize] = 1.0;
        }
    }
    return bands_vec;
}



fn main() {

    //let seq = slim_quick::data_structures::Sequence::new(1234, 4);

    let in_file_name = "/Users/mahdi/Desktop/100.fa".to_string();
    //let in_file_name = "/Users/mahdi/Desktop/seqs.fa".to_string();

    let nb_seqs = 100;
    //let nb_seqs = 1510600;
    let kmer_size = 3;
    let subset_size = 13;
    let nb_subsets = 256usize;
    let nb_kmers = 4usize.pow(kmer_size as u32);


    let max_number_sequences_in_sig = 400;   // max size of a signature before its thrown out

    // Signatures with more max_number_sequences_in_sig sequences are dropped at the end of each iteration
    let mut large_sigs: HashSet<String> = HashSet::new();


    let null_signature = "0.0".to_string();  // handle this better in case it's not 0.0

    let mut signatures = slim_quick::data_structures::Signatures::new();

    let start_time = PreciseTime::now();

    //*************************************************
    // Load sequences in seqs
    //*************************************************


    // load and parse sequences from file
    let mut seqs = data_structures::Sequences::new(nb_seqs);
    seqs.load_sequences_from(in_file_name, kmer_size, nb_subsets);
    let mut curr_time = PreciseTime::now();
    println!("parsing took {} seconds.", start_time.to(curr_time));

    // println!("The loaded sequences are {:?}", seqs);


    //*************************************************
    // Generate sequence counts Matrix
    //*************************************************

    let all_counts = seqs.get_all_counts();
    let counts_matrix = Array::from_shape_vec((nb_seqs, 64), all_counts).unwrap();

    println!("Counts matrix shape size is {:?}", counts_matrix.shape());


    // Add a first column of ONEs that will be multiplied by log_primes, the unique column signature offset
    let counts_augmented_matrix = stack![Axis(1), Array::ones((nb_seqs,1)) , counts_matrix];
    //println!("counts_matrix {:?}", counts_matrix);
    println!("Counts augmented shape size is {:?}", counts_augmented_matrix.shape());

    //*************************************************
    // Generate signatures Matrix
    //*************************************************


    // TODO: clean up and remove this?
    // It was used pre-vectorization
    let bands = utils::get_random_samples(kmer_size, subset_size, nb_subsets, None);

    // Bands_vec is a 0,1 vector of bands of size (nb_kmers x nb_subsets)
    let bands_vec = get_random_samples(nb_kmers, subset_size, nb_subsets, None);
    //*// println!("bands_vec is {:?}", bands_vec);

    // binary marix that will be used for filtering
    // Every column is a mask with exactly subset_size ones and (nb_kemers - subset_size) 0
    let filters_mtx = Array::from_shape_vec((nb_kmers, nb_subsets), bands_vec).unwrap();

    // Column vector of primes
    // we make it of size nb_subset because will need the value for creating column
    // Specific signatures

    // We need this as a column verctor so that we can use it to update the filters matrix
    let logPrimesVec = Array::from_shape_vec((nb_subsets, 1), logPrimesVect(nb_subsets)).unwrap();
    //println!("logPrimesVec\n{:?}", logPrimesVec);


    //*//println!("-- {:?}", logPrimesVec.shape());
    //*//println!("-- {:?}", filters_mtx.shape());

    // Multiple the logPrimes by the 0,1 values
    // we multiply (using boradcasting) to get the log values
    // assigned to the ONEs in the matrix
    let filters_mtx = filters_mtx * &logPrimesVec.slice(s![0..nb_kmers, ..]);
    println!("filters_mtx shape is \n{:?}", filters_mtx.shape());

    // We add a row of logs to guarantee unique signatures
    // This is not going to work... we need logs that haven't been used in the matrix
    // Otherwise we can end up with same sig for seq with diff counts
    let filters_mtx_augmented = stack![Axis(0),  logPrimesVec.slice(s![0..nb_subsets, ..]).t(), filters_mtx];

    println!("filters_mtx augumented shape is \n{:?}", filters_mtx_augmented.shape());


    //*************************************************
    // Do Matrix Multiplication
    //*************************************************


    let start_time = PreciseTime::now();
    let clusters = counts_augmented_matrix.dot(&filters_mtx_augmented);

    //let clusters = clusters.t();

    let mut curr_time = PreciseTime::now();
    println!("matrix multiplication took {} seconds.", start_time.to(curr_time));


    println!("Clusters shape is \n{:?}", clusters.shape());

    //*************************************************
    // parsing resulting clusters matrix
    //*************************************************


    let start_time = PreciseTime::now();

    let mut signatures = data_structures::Signatures::new();
    let mut seen_cluster_definition: HashMultiSet<String> = HashMultiSet::new();
    let mut sig_string_rep: HashSet<String> = HashSet::new();


    let mut clusters_seen = 0;
    // Loop over the clusters, one subset at a time
    for sigs in clusters.axis_iter(Axis(1)) {
        //TOD
        let mut temp_signatures_hash_map: HashMap<String, Vec<i32>> = HashMap::new();
        let mut j = 0;
        for sig in sigs {
            temp_signatures_hash_map.entry(sig.to_string()).or_insert(vec![]).push(j);
            j += 1;
        }
        // Dropping singletons and large signatures
        //println!("{:?}", temp_signatures_hash_map);
        for (sig, seqs_vector) in temp_signatures_hash_map {
            let clusterSize = seqs_vector.len();
            if (1 < clusterSize) && (clusterSize < max_number_sequences_in_sig) {

                // getting string rep of the cluster so that we can
                // test if it already exists
                let seq_strings: Vec<String> = seqs_vector.iter().map(| x | x.to_string()).collect();
                let string_rep = seq_strings.join(",");


                // if signature not seen, then create it.
                if !sig_string_rep.contains(&string_rep) {

                    let sigObj = slim_quick::data_structures::Signature::from_seq_list_with_str_rep(
                        sig.clone(),
                        0,
                        seqs_vector,
                        string_rep);

                    seen_cluster_definition.insert(sigObj.string_rep.clone());
                    signatures.signatures.insert(sig.clone(), sigObj);
                    let sig_str = signatures.signatures[&sig].string_rep.clone();

                    sig_string_rep.insert(sig_str);

                }
            }
        }
        clusters_seen +=1;
        println!("{}", clusters_seen);

    }

    let mut curr_time = PreciseTime::now();
    println!("Parsing of clusters matrix took {} seconds.", start_time.to(curr_time));

    //*************************************************
    // Creating Sequences from Signatures
    //*************************************************
    println!("Started Parsing Signatures");
    let start_time = PreciseTime::now();

    //    for sigObj in signatures.signatures.values() {
    //       println!("12345: {}", sigObj.string_rep);
    //    }

    for (sig, sig_obj) in &signatures.signatures{

        for seq_id in & sig_obj.sequences {
            seqs.sequences[*seq_id as usize].signatures.push(sig.clone());

        }
    }

    let mut curr_time = PreciseTime::now();
    println!("Parsing of signatures took {} seconds.", start_time.to(curr_time));

    //*************************************************
    // Finding the partitons
    //*************************************************

    println!("Started Finding partitions");

    let start_time = PreciseTime::now();

    let my_partitions = slim_quick::find_partitions_v2(nb_seqs, nb_subsets, seqs, &mut signatures);

    let mut curr_time = PreciseTime::now();

    println!("Finding partitions took {} seconds.", start_time.to(curr_time));

    //*************************************************
    // End here
    //*************************************************

    panic!("Temporarily stop here");



    // Where we will temporarily store signatures before processing them
    // to remove bad ones


    let mut subset_number = 0;

    // Equivalent to a Counter with {k -> count for k}
    let mut seen_cluster_definition: HashMultiSet<String> = HashMultiSet::new();

    // we define a cluster as a set of sequences that belong to the same signature
    // ex{seq1, seq72, seq34}
    // If we've already seen a cluster, we will not need to keep track of it again, i.e.,
    // we can discard its signature.
    // We will keep track of a cluster as its sequence representation
    // TODO: Use unique log encoding function here to encode the representation of clusters
    // needs to make sure the sequences in a cluster are sorted.
    let mut string_rep_to_bands: HashMap<String, Vec<u16>> = HashMap::new();
    //let mut string_rep_to_bands_bloomfilter =  Bloom::new(1000000000, 10000000);


    let mut signatures_to_remove: Vec<String> ;



    for subset in bands {

        let mut it_time = PreciseTime::now();

        println!("working on subset {}", subset_number);

        // will be used to remove large stage signatures before we remove
        // large signatures and singletons
        let mut temp_signatures_hash_map:HashMap<String, Vec<i32>> = HashMap::new();
        //        keep track of temp Signatures temp_signatures as
        //            signature => list of sequences
        //        We don't need the subset since we have access to it as a variable in the list
        // The final signatures in temp_signatures_hash_map that make it can then be added the
        // to the permanent signatures list

        let mut signatures_to_remove: Vec<String> = Vec::new();


        let mut nb_sequences_processed = 0;

        for seq in &mut seqs.sequences {

            let sig_value = seq.compute_signature(&subset, subset_number);
            // ignore signature if it contains only empty values
            if sig_value == null_signature {
                continue;
            }

            if !large_sigs.contains(&sig_value) {
                temp_signatures_hash_map.entry(sig_value.clone()).or_insert(vec![]).push(seq.id);

            }else{
                continue;
            }


            if  temp_signatures_hash_map.get(&sig_value).unwrap().len() >= max_number_sequences_in_sig {

                temp_signatures_hash_map.remove(&sig_value);
                large_sigs.insert(sig_value.clone());

                // TODO: ??????
                // continue

            }

            nb_sequences_processed += 1;

            if nb_sequences_processed % 100000 == 0 {
                println!("Processed {} sequences", nb_sequences_processed);
            }

        }// Done processing all the sequencs for this band


        println!("\t\tIn iteration {}, found {} signatures",
                     subset_number, temp_signatures_hash_map.len());
        ////// qprintln!("{:?}", temp_signatures_hash_map);


        // remove signatures containing a single sequence
        let nb_total_sigs_total = temp_signatures_hash_map.len();

        temp_signatures_hash_map.retain(|_, ref v| v.len() > 1 );
        // The number of singletons is the difference in size between old size and new size
        let nb_singletons = nb_total_sigs_total - temp_signatures_hash_map.len();
        println!("\t\tIn iteration {}, removing {} singletons", subset_number, nb_singletons);


        // Removing signatures that contain already seen clusters.
        // TODO: find a better way to do this rather than rely solely on cluster String encoding

//        println!("string_rep_to_bands contains {:?}", string_rep_to_bands_bloomfilter);


        for (sig_id, seqs_vector) in & temp_signatures_hash_map {
            //println!("\t\t\tChecking if Signature was seen: {}", &temp_signature.string_rep);
            let string_seqs_vector:Vec<String> = seqs_vector.iter().map(| &x | x.to_string()).collect();
            let sig_string_rep = string_seqs_vector.join(",");

            if string_rep_to_bands.contains_key(&sig_string_rep)  {
            //if string_rep_to_bands_bloomfilter.check(&temp_signature.string_rep)  {

                // Since the signature was already seen, we can safely remove it here
                signatures_to_remove.push((*sig_id).clone());
                //println!("\t\t\t\tSeen:");

                string_rep_to_bands.get_mut(&sig_string_rep).unwrap().push(subset_number);

            } else {
                string_rep_to_bands.insert(sig_string_rep.clone(), vec![subset_number]);
                //string_rep_to_bands_bloomfilter.set(&temp_signature.string_rep.clone());
                //println!("\t\t\t\tNot seen:");

            }

            seen_cluster_definition.insert(sig_string_rep);
        }


        // Handle all the  signatures_to_remove: remove signature
        // assign the value -1 to all the sequences that have a sequence that will be removed
        println!("\t\tIn iteration {}, removing {} large or seen signatures",
                 subset_number, signatures_to_remove.len());

        for sig_value in &signatures_to_remove{
            temp_signatures_hash_map.remove(sig_value);
        }

        println!("The Number of temp signatures is {}", temp_signatures_hash_map.len());
        // println!("\n\n\t\tseen_cluster_definition is: {:?}", seen_cluster_definition.distinct_elements());
        // println!("\t\tstring_rep_to_bands is: {:?}", string_rep_to_bands);


        println!("--------------------------");

        // Here all the temp_signatures_hash_map seem legit. Create and add them to signatures collection;
        for (sig_id, seqs_vector) in temp_signatures_hash_map {

            let sig = slim_quick::data_structures::Signature::from_seq_list(sig_id.clone(), 0, seqs_vector);

            signatures.signatures.insert(sig_id, sig);


        }


//        for (sig_id, signature) in &signatures.signatures{
//            for tempSeq in &signature.sequences {
//                 seqs.sequences[*tempSeq as usize].signatures[signature.subset_number as usize] = sig_id.clone();
//            }
//        }


        // signatures.extend(temp_signatures);
        println!("\t\tTotal number of signatures so far is {}", signatures.len());

        let mut time_now = PreciseTime::now();
        println!("iteration {} took {} seconds.", subset_number, it_time.to(time_now));

        subset_number += 1;
        // temp_signatures.clear() was automatically



    }


    // make sure, at least for testing purposes, that the number of elems in
    // string_rep_to_bands is the same as the number of elements in seen_cluster_definition
    // TODO: Remove this later or run only when testing
    //assert_eq!(string_rep_to_bands.len(), seen_cluster_definition.distinct_elements().collect::<HashSet<_>>().len());

    // println!("After all iterations, total number of clusters is {}", string_rep_to_bands.len());

    // return string_reps that were only seen once.
    // string_rep_to_bands = string_rep_to_bands.into_iter()
    //      .filter(| &(_, ref y) |  y.len() == 1)
    //      .map(|(k,v)| (k,v))
    //      .collect();

    // TODO: remove clusters that were only seen once?
    //println!("After removing all clusters only seen once, total number of clusters is {}", string_rep_to_bands.len());

    let mut partition_time = PreciseTime::now();

    for x in seen_cluster_definition.distinct_elements(){
        println!("{}::{}", x, seen_cluster_definition.count_of(x.clone()));

    }
    println!("\n\n\n----------------------------------------\n\n\n");
    println!("Before removing bad signatures {:?}", seqs.sequences);
    for (sigString, sig) in &signatures.signatures{

        if seen_cluster_definition.count_of(sig.string_rep.clone()) > 1{

            println!("{}", sig.string_rep);
            println!("{}", seen_cluster_definition.count_of(sig.string_rep.clone()));

            for tempSeq in &signatures.signatures[sigString].sequences {
                seqs.sequences[*tempSeq as usize].signatures.push(sigString.clone());
            }

        }
    }

    println!("After removing bad signatures {:?}", seqs.sequences);

    println!("my Seen clusters are");

    println!("Finding the partitons");
    //let my_partitions = slim_quick::find_partitions_v2(nb_seqs, nb_subsets, seqs, &mut signatures);
//    // println!("My partitions are: {:?}", my_partitions);
//
//    let mut time_now = PreciseTime::now();
//    println!("Finding partition took {} seconds.", partition_time.to(time_now));
//
//
//
//    curr_time  = PreciseTime::now();
//
//    println!("Total runtime is {} seconds.", start_time.to(curr_time));



}
