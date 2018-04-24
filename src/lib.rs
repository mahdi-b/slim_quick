#![allow(non_snake_case)]
#![allow(unused_variables)]
#![allow(dead_code)]

#[macro_use]



extern crate lazy_static;
extern crate rand;
extern crate bio;

use std::collections::BinaryHeap;
use std::collections::HashSet;





pub mod data_structures;
pub mod utils;



pub fn find_partitions(nb_sequences: usize, nb_bands: usize,
                       seqs: data_structures::Sequences, sigs: &mut data_structures::Signatures)
                       -> Vec<i32> {
    //TODO check that len of seqs == nb_sequences

    println!("I have passed along to find_partitions {} sequences and {} signatures", seqs.sequences.len(), sigs.signatures.len());

    //println!("{:?}", sigs.signatures);

    //println!("{:?}", seqs);

    let mut partitions: Vec<i32>  = vec![-1; nb_sequences as usize];

    let mut my_queue: BinaryHeap<usize> =  BinaryHeap::new();

    let mut neighbors_set: HashSet<i32> = HashSet::new();
    //let mut neighborhood_count: HashMap<u32, u16> = HashMap::new();

    let mut next_partiton_id = 0; // nbPartitions > max int value??

    for i  in 0..nb_sequences {

        if partitions[i] == -1 {
            // Assign it to its own cluster
            // Cannot belong to a previous cluster, or else it would have been labeled through another sequence
            partitions[i] = next_partiton_id;
            next_partiton_id+=1;
        } else {
            // Was already assigned as part of another sequence
            continue;
        }
        // Add the sequence to a queue so as to inspect its children
        my_queue.push(i);

        while ! my_queue.is_empty() {

            let new_seq_id:i32 = my_queue.pop().unwrap() as i32;

            // For each signature that this sequence is involved in
            // get all the other sequences that are in the same signatures
            for band_number in 0..nb_bands {

                let sig_id:String = seqs.sequences[new_seq_id as usize].signatures[band_number].to_string();

//                println! ("Signatures for seq_id {} are {:?}", new_seq_id, seqs.sequences[new_seq_id as usize].signatures);
//                println!("Inspecting sig_id {} ", sig_id);
//                println!("Identifying the neighbors of {:?} ", sig.sequences);

                if sig_id == "-1"{
                    // signature was deleted of has already been visited
                    continue;
                }
                let sig = sigs.signatures.get_mut(&sig_id).unwrap();



                if sig.visited == true {
                    // signature was deleted of has already been visited
                    continue;
                }




                // TODO: fix to start at the index of the sequence
                // We start at the position where we found the sequence,  since all sequences
                // before it would have been assigned already.
                // Given that seqIds are incremental we are sure
                // that any sequence before it has already been assigned to a partition
                for seq_idx in 0..sig.sequences.len() {
                    let seq_sig = sig.sequences[seq_idx];
                    if seq_sig != new_seq_id && partitions[seq_sig as usize] == -1 {
                        // Sequence was not assigned to a cluster
                        // so add it to the queue to process other sequences in its cluster
                        neighbors_set.insert(seq_sig);
                    }
                }
                sig.visited = true;
            }
            // Now that we have all the neighbors in the same signature as this one,
            // Assign them all to the same signature as sequence i
            // and add them to the queue so that we can process their neighbors
            for neighbor_of_i in neighbors_set.iter() {
                partitions[*neighbor_of_i as usize] = partitions[i as usize];
                // cout << *it << endl;
                my_queue.push(*neighbor_of_i as usize);
            }

            neighbors_set.clear();

        }
    }
    println!("partitions is {:?}", partitions);

    return partitions;
}



// This version finds partitions based on the temp_signatures_hash_map
pub fn find_partitions_v2(nb_sequences: usize, nb_bands: usize,
                       seqs: data_structures::Sequences, sigs: &mut data_structures::Signatures)
                       -> Vec<i32> {
    //TODO check that len of seqs == nb_sequences

    //println!("I have passed along to find_partitions {} sequences and {} signatures", seqs.sequences.len(), sigs.signatures.len());

    //println!("{:?}", sigs.signatures);

    //println!("{:?}", seqs);

    let mut partitions: Vec<i32>  = vec![-1; nb_sequences as usize];

    let mut my_queue: BinaryHeap<usize> =  BinaryHeap::new();

    let mut neighbors_set: HashSet<i32> = HashSet::new();
    //let mut neighborhood_count: HashMap<u32, u16> = HashMap::new();

    let mut next_partiton_id = 0; // nbPartitions > max int value??

    for i  in 0..nb_sequences {

        if partitions[i] == -1 {
            // Assign it to its own cluster
            // Cannot belong to a previous cluster, or else it would have been labeled through another sequence
            partitions[i] = next_partiton_id;
            next_partiton_id+=1;
        } else {
            // Was already assigned as part of another sequence
            continue;
        }
        // Add the sequence to a queue so as to inspect its children
        my_queue.push(i);

        while ! my_queue.is_empty() {

            let new_seq_id:i32 = my_queue.pop().unwrap() as i32;

            // For each signature that this sequence is involved in
            // get all the other sequences that are in the same signatures
            for sig_id in &seqs.sequences[new_seq_id as usize].signatures {


//                println! ("Signatures for seq_id {} are {:?}", new_seq_id, seqs.sequences[new_seq_id as usize].signatures);
//                println!("Inspecting sig_id {} ", sig_id);
//                println!("Identifying the neighbors of {:?} ", sig.sequences);

                if sig_id == "-1"{
                    // signature was deleted of has already been visited
                    continue;
                }
                let sig = sigs.signatures.get_mut(sig_id).unwrap();



                if sig.visited == true {
                    // signature was deleted of has already been visited
                    continue;
                }




                // TODO: fix to start at the index of the sequence
                // We start at the position where we found the sequence,  since all sequences
                // before it would have been assigned already.
                // Given that seqIds are incremental we are sure
                // that any sequence before it has already been assigned to a partition
                for seq_idx in 0..sig.sequences.len() {
                    let seq_sig = sig.sequences[seq_idx];
                    if seq_sig != new_seq_id && partitions[seq_sig as usize] == -1 {
                        // Sequence was not assigned to a cluster
                        // so add it to the queue to process other sequences in its cluster
                        neighbors_set.insert(seq_sig);
                    }
                }
                sig.visited = true;
            }
            // Now that we have all the neighbors in the same signature as this one,
            // Assign them all to the same signature as sequence i
            // and add them to the queue so that we can process their neighbors
            for neighbor_of_i in neighbors_set.iter() {
                partitions[*neighbor_of_i as usize] = partitions[i as usize];
                // cout << *it << endl;
                my_queue.push(*neighbor_of_i as usize);
            }

            neighbors_set.clear();

        }
    }
    println!("partitions is {:?}", partitions);

    return partitions;
}



mod test {
    mod test_find_partitions {
        use super::super::*;
        #[test]
        fn test_find_partitions(){
            let mut seq1: data_structures::Sequence = data_structures::Sequence::new(1, vec![3,0,2,3,1], 3);
            let mut seq2: data_structures::Sequence = data_structures::Sequence::new(2, vec![3,0,2,3,1], 3);
            let mut seq3: data_structures::Sequence  = data_structures::Sequence::new(3, vec![4,0,2,1,0], 3);

            seq1.add_signature_for_band("1".to_string(), 0);
            seq1.add_signature_for_band("2".to_string(), 1);
            seq1.add_signature_for_band("3".to_string(), 2);

            seq2.add_signature_for_band("1".to_string(), 0);
            seq2.add_signature_for_band("2".to_string(), 1);
            seq2.add_signature_for_band("3".to_string(), 2);

            seq3.add_signature_for_band("4".to_string(), 0);
            seq3.add_signature_for_band("5".to_string(), 1);
            seq3.add_signature_for_band("6".to_string(), 2);

            let mut seqs = data_structures::Sequences::new(3);
            seqs.add_seq(seq1);
            seqs.add_seq(seq2);
            seqs.add_seq(seq3);

            let mut sig1 = data_structures::Signature::new("1".to_string(), 500, 0);
            let mut sig2 = data_structures::Signature::new("2".to_string(), 500, 1);
            let mut sig3 = data_structures::Signature::new("3".to_string(), 500, 2);
            let mut sig4 = data_structures::Signature::new("4".to_string(), 500, 0);
            let mut sig5 = data_structures::Signature::new("5".to_string(), 500, 1);
            let mut sig6 = data_structures::Signature::new("6".to_string(), 500, 2);

            sig1.add_sequence(0);
            sig1.add_sequence(1);

            sig2.add_sequence(0);
            sig2.add_sequence(1);

            sig3.add_sequence(0);
            sig3.add_sequence(1);
            //sig3.add_sequence(2);


            sig4.add_sequence(2);

            sig5.add_sequence(2);

            sig6.add_sequence(2);


            let mut sigs = data_structures::Signatures::new();
            sigs.signatures.insert("1".to_string(), sig1);
            sigs.signatures.insert("2".to_string(), sig2);
            sigs.signatures.insert("3".to_string(), sig3);
            sigs.signatures.insert("4".to_string(), sig4);
            sigs.signatures.insert("5".to_string(), sig5);
            sigs.signatures.insert("6".to_string(), sig6);




            find_partitions_v2(3, 3, seqs, &mut sigs);


//            assert_eq!(1,0);


        }

    }
}
