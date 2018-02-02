#![allow(unused_variables)]
#![allow(unused_mut)]
#![allow(unused_imports)]
//TODO: replace the word band with subset throughout the code

extern crate slim_quick;

//mod data_structures;
//mod utils;
use slim_quick::data_structures;
use slim_quick::utils;


#[macro_use]
extern crate lazy_static;

extern crate rand;

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





fn main(){

    //let seq = slim_quick::data_structures::Sequence::new(1234, 4);

    let in_file_name = "/Users/mahdi/Desktop/seqs.fna".to_string();

    let nb_seqs = 1510600;
    let kmer_size = 3;
    let subset_size = 15;
    let nb_subsets = 100usize;


    let max_number_sequences_in_sig = 500;   // max size of a signature before its thrown out

    // Signatures with more max_sig_size sequences are dropped at the end of each iteration
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

    // to keep track of the bands in which a cluster was seen
    // if we've already seen a cluster, we will not need to keep track of it again
    // we will keep track of a cluster as its sequence representation
    // TODO: Use unique log encoding function here to encode the representation of clusters
    // needs to make sure the sequences in a cluster are sorted.
    let mut string_rep_to_bands: HashMap<String, Vec<u16>> = HashMap::new();

    let mut signatures_to_remove: Vec<String> ;



    for subset in bands {

        let mut it_time = PreciseTime::now();

        println!("working on subset {}", subset_number);

        // will be used to remove large stage signatures before we remove
        // large signatures and singletons
        let mut temp_signatures = data_structures::Signatures::new();
        let mut signatures_to_remove: Vec<String> = Vec::new();


        let mut nb_sequences_processed = 0;

        for seq in &mut seqs.sequences {
            let sig_value = seq.compute_signature(&subset, subset_number);


            // ignore signature if it contains only empty values
            if sig_value == null_signature {
                continue;
            }


            if !large_sigs.contains(&sig_value) {
                // create the signature so that we can add it.
                // better if we create the signature in or_insert()


                if ! temp_signatures.signatures.contains_key(&sig_value){
                    let mut signature: data_structures::Signature =
                        data_structures::Signature::new(sig_value.clone(),
                                                        max_number_sequences_in_sig,
                                                        subset_number);
                    temp_signatures.signatures.insert(sig_value.clone(), signature);
                }
            } else {
                // signature is in large_sigs already so
                // this sequence won't have a signature for this subset
                // seq.sigIds[subsetNumber] = -1 as initialized
                continue;
            }


            let signature = temp_signatures.signatures.get_mut(&sig_value).unwrap();

            // if signature not too large
            if signature.sequences.len() < max_number_sequences_in_sig {

                signature.add_sequence(seq.id);
                // TODO: instead of cloning signature_id, we can store a ref to the signature id
                seq.add_signature_for_band(signature.signature_id.clone(), subset_number);

                // String rep(resentation) is a comma del list of sequences_id
                if signature.string_rep.len() > 0 {
                    signature.string_rep.push_str(&format!(", {}", seq.id));
                }else{
                    signature.string_rep.push_str(&format!("{}", seq.id));
                }
            } else {
                // remove the signature but also keep track of it so that we don't add to
                // it again
                signatures_to_remove.push(signature.signature_id.clone());

                // if !large_sigs.contains(&signature.signature_id) {
                large_sigs.insert(signature.signature_id.clone());

            }

            nb_sequences_processed += 1;

            if nb_sequences_processed % 100000 == 0 {
                println!("Processed {} sequences", nb_sequences_processed);
            }

        }// Done processing all the sequencs for this band


        // Just for status update purposes
        if temp_signatures.signatures.len() < 10 {
            println!("\t\tIn iteration {}, found {} signatures \n\t\t {:?}",
                     subset_number, temp_signatures.signatures.len(), temp_signatures.signatures);
        }else{
            println!("\t\tIn iteration {}, found {} signatures",
                     subset_number, temp_signatures.signatures.len());

        }



        // remove signatures containing a single sequence (singletons signatures) and sets value of sequences
        //containing that signature to -1
        let singleton_signatures = temp_signatures.remove_singletons();

        println!("\t\tIn iteration {}, removing {} singletons", subset_number, singleton_signatures.len());
        for signature in singleton_signatures {
            seqs.remove_signature_from_sequences(signature)
        }


        // Removing signatures that contain already seen clusters.
        // TODO: find a better way to do this rather tha rely solely on cluster String encoding.


        for (sig_string , temp_signature) in &temp_signatures.signatures {

            if string_rep_to_bands.contains_key(&temp_signature.string_rep)  {

                // Since the signature was already seen, we can safely remove it here
                signatures_to_remove.push(temp_signature.signature_id.clone());


                string_rep_to_bands.get_mut(&temp_signature.string_rep).unwrap().push(subset_number);

            } else {
                string_rep_to_bands.insert(temp_signature.string_rep.clone(), vec![subset_number]);
            }
            seen_cluster_definition.insert(temp_signature.string_rep.clone());
        }


        // Handle all the  signatures_to_remove: remove signature
        // assign the value -1 to all the sequences that have a sequence that will be removed
        println!("\t\tIn iteration {}, removing {} large or seen signatures",
                 subset_number, signatures_to_remove.len());
        for sig_value in &signatures_to_remove{
            let signature: data_structures::Signature = temp_signatures.signatures.remove(sig_value).unwrap();
            seqs.remove_signature_from_sequences(signature);
        }


        // println!("\n\n\t\tseen_cluster_definition is: {:?}", seen_cluster_definition.distinct_elements());

        // println!("\t\tstring_rep_to_bands is: {:?}", string_rep_to_bands);

        println!("--------------------------");

        signatures.extend(temp_signatures);
        println!("\t\tTotal number of signatures so far is {}", signatures.len());

        let mut time_now = PreciseTime::now();
        println!("iteration {} took {} seconds.", subset_number, it_time.to(time_now));


        subset_number += 1;
        // temp_signatures.clear() was automatically

    }



    // make sure, at least for testing purposes, that the number of elems in
    // string_rep_to_bands is the same as the number of elements in seen_cluster_definition
    // TODO: Remove this later or run only when testing
    assert_eq!(string_rep_to_bands.len(), seen_cluster_definition.distinct_elements().collect::<HashSet<_>>().len());

    println!("After all iterations, total number of clusters is {}", string_rep_to_bands.len());

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
    // println!("My partitions are: {:?}", my_partitions);

    let mut time_now = PreciseTime::now();
    println!("Finding partition took {} seconds.", partition_time.to(time_now));



    curr_time  = PreciseTime::now();

    println!("Total runtime is {} seconds.", start_time.to(curr_time));



}
