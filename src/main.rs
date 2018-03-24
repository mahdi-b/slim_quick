#![allow(unused_variables)]
#![allow(unused_mut)]
#![allow(unused_imports)]
//TODO: replace the word band with subset throughout the code

extern crate slim_quick;


extern crate ndarray;

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



fn main(){

    //let seq = slim_quick::data_structures::Sequence::new(1234, 4);

    let in_file_name = "/Users/mahdi/Desktop/100.fa".to_string();

    let nb_seqs = 100;
    let kmer_size = 3;
    let subset_size = 15;
    let nb_subsets = 100usize;


    let max_number_sequences_in_sig = 500;   // max size of a signature before its thrown out

    // Signatures with more max_number_sequences_in_sig sequences are dropped at the end of each iteration
    let mut large_sigs: HashSet<String> = HashSet::new();


    let null_signature = "0.0".to_string();  // handle this better in case it's not 0.0

    let mut signatures  = slim_quick::data_structures::Signatures::new();

    let start_time = PreciseTime::now();


    // load and parse sequences from file
    let mut seqs = data_structures::Sequences::new(nb_seqs);
    seqs.load_sequences_from(in_file_name, kmer_size);
    let mut curr_time = PreciseTime::now();
    println!("parsing took {} seconds.", start_time.to(curr_time));




    // generate bands (subsets)
    let bands = utils::get_random_samples(kmer_size, subset_size, nb_subsets, None);
    println!("bands are {:?}", bands);
    //panic!("I am prepamuruly dying to print output");

    // Where we will temporarily store signatures before processing them
    // to remove bad ones
    let mut signatures = data_structures::Signatures::new();



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
        println!("{:?}", temp_signatures_hash_map);


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

            seen_cluster_definition.insert(sig_string_rep.clone());
        }


        // Handle all the  signatures_to_remove: remove signature
        // assign the value -1 to all the sequences that have a sequence that will be removed
        println!("\t\tIn iteration {}, removing {} large or seen signatures",
                 subset_number, signatures_to_remove.len());

        for sig_value in &signatures_to_remove{
            temp_signatures_hash_map.remove(sig_value);
        }

        println!("2 Number of temp signatures is {}", temp_signatures_hash_map.len());
        // println!("\n\n\t\tseen_cluster_definition is: {:?}", seen_cluster_definition.distinct_elements());
        // println!("\t\tstring_rep_to_bands is: {:?}", string_rep_to_bands);


        println!("--------------------------");
        // Here all the temp_signatures_hash_map seem legit. Create and add them to signatures collection;
        for (sig_id, seqs_vector) in temp_signatures_hash_map {
           let sig = slim_quick::data_structures::Signature::from_seq_list(sig_id.clone(), 4, seqs_vector);
            signatures.signatures.insert(sig_id, sig);
        }

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

    println!("Finding the partitons");
    let my_partitions = slim_quick::find_partitions(nb_seqs, nb_subsets, seqs, &mut signatures);
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
