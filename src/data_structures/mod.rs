use std::process;
use bio::io::fasta;
use bio::io::fasta::Record;
use std::collections::HashMap;

use utils;


lazy_static! {
    pub static ref UTILS_INSTANCE: utils::Utils = utils::Utils::new();
}

//////////////////////////////////////////////////////
//  MARK Struct Signature
//////////////////////////////////////////////////////
#[derive(Clone, Debug)]
pub struct Signature{
    pub signature_id: String,
    pub visited: bool,
    pub string_rep: String, // string representation from sequences in this signature
    pub sequences: Vec<i32>,
    pub subset_number:u16 // subset number to which this signature belongs
}

impl Signature{

    pub fn new(signature_id: String, max_num_sequences:usize, subset_number: u16) -> Signature {
        Signature{
            signature_id: signature_id,
            visited: false,
            string_rep: "".to_string(),
            sequences: Vec::with_capacity(max_num_sequences),
            subset_number:subset_number
        }
    }

    pub fn from_seq_list(signature_id: String, subset_number: u16, sequences: Vec<i32>) -> Signature {

        let seq_strings: Vec<String> = sequences.iter().map(| x | x.to_string()).collect();
        let string_rep:String = seq_strings.join(",");

        Signature {
            signature_id: signature_id,
            visited: false,
            string_rep: string_rep,
            sequences: sequences,
            subset_number: subset_number
        }
    }


    pub fn add_sequence(&mut self, seq_id:i32){
        self.sequences.push(seq_id);
        //println!("added sequence {} to signature {}", seq_id, self.signature_id);
    }

    pub fn get_pos_seq_id(&self, seq_id: i32) -> Option<usize> {
        self.sequences.binary_search(&seq_id).ok()
    }
}


pub struct Signatures{
    pub signatures: HashMap<String, Signature>
}

impl Signatures{
    pub fn new() -> Signatures {
        let sigs:  HashMap<String, Signature> = HashMap::new();
        Signatures{signatures:sigs}
    }


    // TODO return an option instead
    // the returns Vec contains are singletons
    pub fn remove_singletons(&mut self) -> Vec<Signature> {
        let singletons: Vec<_> = self.signatures.iter()
            .filter(|&(k,v)| v.sequences.len() == 1)
            .map(|(k, v)| v.clone())
            .collect();
        //singletons.abra();
        for singleton in  &singletons { self.signatures.remove(&singleton.signature_id); }
        singletons
    }

    pub fn extend(&mut self, temp_signatures: Signatures){
        self.signatures.extend(temp_signatures.signatures)
    }




    pub fn clear(&mut self){
        self.signatures.clear();
    }
    pub fn len(&mut self) -> usize{
        self.signatures.len()
    }


}


//////////////////////////////////////////////////////
// MARK Struct kmer
//////////////////////////////////////////////////////

pub struct Kmer {

}

impl Kmer{

    // getHash for DNA kmer
    pub fn get_encoding(kmer: &[u8]) -> usize {
        let mut hash_code:usize = 0;
        for base in kmer {
            hash_code <<= 2;
            if let Some(val) = utils::convert_byte_DNA(base) {
                hash_code |= val;
            }else{
                println!("Something went wrong!");
            }
        }
        return hash_code;
    }
}


//////////////////////////////////////////////////////
// MARK Struct kmerIter
//////////////////////////////////////////////////////
pub struct KmerIter{
    position: usize,
    kmer_size: usize,
    sequence_buffer: Vec<u8>
}

impl Iterator for KmerIter{
    type Item = (usize,  String);

    fn next(&mut self) -> Option<(usize, String)> {

        if (self.position  + self.kmer_size) > self.sequence_buffer.len() {
            return None
        }
        self.position += 1;
        let kmer_vec =  self.sequence_buffer[self.position-1 .. self.position-1+self.kmer_size].to_vec();

        let kmer = String::from_utf8(kmer_vec).expect("found unknown char in kmer encoding");
        Some((self.position-1, kmer))
    }
}


//////////////////////////////////////////////////////
// MARK Struct Sequence
//////////////////////////////////////////////////////
#[derive(Debug)]
pub struct Sequence{
    pub id: i32,
    pub kmers_counts: Vec<u16>,
    pub signatures: Vec<String>
}

impl Sequence{


    pub fn new (seq_id: i32, kmers_counts: Vec<u16>,  nb_signatures: usize) -> Sequence {
        print!("creating a new sequence with seq_id{}", seq_id);
        Sequence{ id: seq_id, kmers_counts: kmers_counts, signatures: Vec::with_capacity(nb_signatures) }
    }

    // making get_kmers static so as not to retain the sequence in memory
    pub fn get_kmers(sequence: String, kmer_size:usize) -> KmerIter{
        let sequence_vec = sequence.as_bytes().to_vec();
        KmerIter{position:0, kmer_size:kmer_size, sequence_buffer:sequence_vec}
    }

    pub fn get_kmers_counts(sequence: String, kmer_size:usize) -> Vec<u16>{
        let mut counts = vec![0u16; 4usize.pow(kmer_size as u32) ];
        for (_, kmer) in Sequence::get_kmers(sequence, kmer_size){
            counts[Kmer::get_encoding(kmer.as_bytes())] +=1;
        }
        counts
    }




    pub fn compute_signature(&self, band: &[u16], subsetNumber:u16) -> String{

        let mut signature: f64=0.0;
        // println!("bands are: {:?}", band);
        for col in band{
            //println!("multiplying {}  and {}", self.kmers_counts[*col as usize]as f64, UTILS_INSTANCE.log_primes[*col as usize]);
            signature +=   (self.kmers_counts[*col as usize] as f64) * UTILS_INSTANCE.log_primes[*col as usize];
        }
        //println!("multiplying {}  and {}", (subsetNumber as f64), UTILS_INSTANCE.log_primes[256]);

        signature +=   (subsetNumber as f64) * UTILS_INSTANCE.log_primes[256];
        signature.to_string()
    }

    pub fn add_signature_for_band(&mut self, signature_id: String, bandNumber: u16){
        //TODO, we are not using band number here.... remove it?
        self.signatures.push(signature_id)
    }

    pub fn get_signature_for_band(&self, bandNumber: u16) -> Option<&String>{
        if (bandNumber as usize) < self.signatures.len() {
            Some(&self.signatures[bandNumber as usize])
        }else{
            None
        }
    }
    // TODO Add tests for this
    pub fn remove_signature_for_band(&mut self, bandNumber: u16) {
        self.signatures[bandNumber as usize] = "-1".to_string();
    }
}


//////////////////////////////////////////////////////
// MARK Struct Sequences
//////////////////////////////////////////////////////
#[derive(Debug)]
pub struct Sequences{
    pub sequences: Vec<Sequence>,


}


impl Sequences {
    pub fn new(nb_sequences: usize) -> Sequences{
        let seqs = Sequences{sequences: Vec::with_capacity(nb_sequences)};
        println!("Creating the list of sequences");
        seqs
    }

    pub fn load_sequences_from(& mut self, in_file_name: String, kmer_size: usize, nb_subsets: usize) {

        let reader = fasta::Reader::from_file(& in_file_name).unwrap_or_else(|err|{
            println!("Problem opening input file {}: {}", in_file_name, err);
            process::exit(1);
        });

        println!("Parsing sequences from file {} ", in_file_name);
        let mut i= 0;
        for record in reader.records(){
            let record: Record = record.unwrap();
            let seq = Sequence {
                id: record.id().unwrap().parse().unwrap(),
                kmers_counts: Sequence::get_kmers_counts(String::from_utf8(record.seq().to_vec()).unwrap(), kmer_size),
                signatures: Vec::new()
                    // signatures: vec!["-1".to_string(); nb_subsets]
            };

            self.add_seq(seq);

            i+=1 ;
            if i % 100000 == 0{
                println!("\t reads {} sequences", i);
            }

            // println!("Sequence: {:?}", record.id());
            // println!("{:?}", record.seq());
        }
        println!("Done parsing sequences from file {} ", in_file_name);
    }

    pub fn add_seq(& mut self, seq: Sequence){
        self.sequences.push(seq);
    }

    pub fn remove_signature_from_sequences(&mut self, signature: Signature){
        // TODO: fix this cloned thing below
        let missing = "-1".to_string();
        for seq_id in signature.sequences{
            self.sequences[seq_id as usize].signatures[signature.subset_number as usize] = missing.clone();
        }
    }

    pub fn get_all_counts(&mut self) -> Vec<f64>{

        let mut all_counts = vec![];

        for seq in &self.sequences{
            for c in &seq.kmers_counts {
                all_counts.push(*c as f64);
            }
        }
        all_counts
    }



}




//////////////////////////////////////////////////////
// Mark Tests
//////////////////////////////////////////////////////

#[cfg(test)]
mod test{

    //////////////////////////////////////////////////////
    // MARK KmerIter tests
    //////////////////////////////////////////////////////
    mod kmerIter_tests {
        //use super::super::*;
    }

    //////////////////////////////////////////////////////
    // MARK signature tests
    //////////////////////////////////////////////////////
    mod signature_tests {
        use super::super::*;

        #[test]
        fn test_new(){
            let signature_1  = "1234".to_string();
            let sig =  Signature::new(signature_1, 4, 0);
            assert_eq!(sig.sequences.len(),0);
            assert_eq!(sig.visited,false);
            assert_eq!(sig.signature_id, "1234".to_string());

        }
        #[test]
        fn test_from_seq_list(){

            let sig_id  = "1234".to_string();
            let sequences = vec![0,1,2,3];
            let sig = Signature::from_seq_list(sig_id, 4, sequences);
            assert_eq!(sig.sequences.len(),4);
            assert_eq!(sig.sequences[0], 0);
            assert_eq!(sig.sequences[3], 3);
            assert_eq!(sig.signature_id, "1234".to_string());


        }


        #[test]
        fn test_add_sequence() {
            let signature_1  = "1234".to_string();
            let mut sig = Signature::new(signature_1, 4, 0);
            sig.add_sequence(1);
            assert_eq!(sig.sequences.len(), 1);
        }
        #[test]
        fn test_get_pos_seq_id() {
            let signature_1  = "1234".to_string();
            let mut sig = Signature::new(signature_1, 4, 0);
            sig.add_sequence(1);
            sig.add_sequence(2);
            sig.add_sequence(3);
            sig.add_sequence(4);
            assert_eq!(sig.get_pos_seq_id(1), Some(0));
            assert_eq!(sig.get_pos_seq_id(4), Some(3));
        }
    }

    //////////////////////////////////////////////////////
    // MARK signature tests
    //////////////////////////////////////////////////////
    mod signatures_tests {

        use super::super::*;

        #[test]
        fn test_new(){
            let signatures = Signatures::new();
            assert_eq!(signatures.signatures.len(), 0)
        }
        #[test]
        fn test_remove_singletons(){

            let mut sigs = Signatures::new();

            let mut  sig =  Signature::new("123".to_string(), 4, 0);
            sig.add_sequence(0);
            sig.add_sequence(1);
            sigs.signatures.insert(sig.signature_id.clone(), sig);


            let mut  sig =  Signature::new("456".to_string(), 4, 0);
            sig.add_sequence(0);
            sig.add_sequence(1);
            sigs.signatures.insert(sig.signature_id.clone(), sig);


            let mut  sig =  Signature::new("789".to_string(), 4, 0);
            sig.add_sequence(0);
            sigs.signatures.insert(sig.signature_id.clone(), sig);

            sigs.remove_singletons();
            assert_eq!(sigs.signatures.len(), 2);
            assert!(! sigs.signatures.contains_key(&"789".to_string()));
        }
    }


    //////////////////////////////////////////////////////
    // MARK kmer tests
    //////////////////////////////////////////////////////
    mod kmer_tests {
        use super::super::*;
        #[test]
        fn test_get_kmer_encoding(){
            assert_eq!(Kmer::get_encoding(b"A"), 0);
            assert_eq!(Kmer::get_encoding(b"C"), 1);
            assert_eq!(Kmer::get_encoding(b"G"), 2);
            assert_eq!(Kmer::get_encoding(b"T"), 3);
            assert_eq!(Kmer::get_encoding(b"AA"), 0);
            assert_eq!(Kmer::get_encoding(b"CA"), 4);
            assert_eq!(Kmer::get_encoding(b"GA"), 8);
            assert_eq!(Kmer::get_encoding(b"TA"), 12);
            assert_eq!(Kmer::get_encoding(b"ACGT"), 27);
            assert_eq!(Kmer::get_encoding(b"TGCA"), 228);
        }
    }



    //////////////////////////////////////////////////////
    // MARK sequence tests
    //////////////////////////////////////////////////////
    mod sequence_tests {
        use super::super::*;
        #[test]
        fn test_get_kmers() {
            for (pos, kmer) in Sequence::get_kmers("ACGTA".to_string(), 3) {
                if pos == 0 {
                    assert_eq!(kmer, "ACG");
                } else if pos == 1 {
                    assert_eq!(kmer, "CGT");
                    //assert!(false);
                }else if pos == 2 {
                    assert_eq!(kmer, "GTA");
                } else {
                    println!("Unknown value for pos {}", pos);
                    assert!(false);
                }
            }
        }

        #[test]
        fn test_get_kmers_counts(){

            let kmer_size = 3;
            let seq = "ACGTACGTCGT".to_string();
            let counts = Sequence::get_kmers_counts(seq, kmer_size);

            assert_eq!(counts[Kmer::get_encoding(b"CGT")] , 3);
            assert_eq!(counts[Kmer::get_encoding(b"ACG")] , 2);
            assert_eq!(counts[Kmer::get_encoding(b"GTA")] , 1);
            assert_eq!(counts[Kmer::get_encoding(b"AAA")] , 0);
            //println!("{:?}", counts.iter().fold(0,|a, &b| a + b));
        }

        #[test]
        fn test_new(){
            let seq = Sequence::new(1, vec![0,1,2,3,4], 4);
            assert_eq!(seq.id, 1);
            assert_eq!(seq.signatures.capacity(), 4);
        }


        #[test]
        fn test_kmers(){
            // TODO Add a test for this
            assert_eq!(1, 1);
        }

        #[test]
        fn test_compute_signature(){
            let seq = Sequence::new(1, vec![3,0,2,3,1], 4);
            let mut signature:f64  = seq.compute_signature(&[0,1,4], 1).parse().unwrap();
            // vals (3,0,1), primes are (0.3010299956639812, 0.47712125471966244, 1.0413926851582251)
            // subset_number 1 and UTILS_INSTANCE.log_primes[256] = 3.209783014848515
            assert!(( signature - 5.1542656869986843).abs() < 0.00001);
            signature  = seq.compute_signature(&[3,4], 5).parse().unwrap();
            // vals (3,1), primes are ( 0.8450980400142568, 1.0413926851582251)
            // subset_number 5 and UTILS_INSTANCE.log_primes[256] = 3.209783014848515
            assert!((signature - 19.625601879443572).abs() < 0.00001);
        }



        #[test]
        fn test_add_signature(){
            let signature_1  = "1234".to_string();
            let signature_2  = "5678".to_string();

            let sig1 = Signature::new(signature_1, 4, 0);
            let sig2 = Signature::new(signature_2, 4, 0);

            let mut seq = Sequence::new(1, vec![0,1,2,3,4], 4);
            seq.add_signature_for_band(sig1.signature_id, 0);
            seq.add_signature_for_band(sig2.signature_id, 1);
            assert_eq!(*seq.signatures[1], "5678".to_string());
        }

        #[test]
        fn test_get_signature(){
            let signature_1  = "1234".to_string();
            let signature_2  = "5678".to_string();

            let sig1 = Signature::new(signature_1, 4, 0);
            let sig2 = Signature::new(signature_2, 4, 0);

            let mut seq = Sequence::new(1, vec![0,1,2,3,4], 4);
            seq.add_signature_for_band(sig1.signature_id, 0);
            seq.add_signature_for_band(sig2.signature_id, 1);

            println!("******************************\n******************************\n");
            assert_eq!(*seq.get_signature_for_band(1).unwrap(), "5678".to_string());
            println!("{:?}", seq.get_signature_for_band(10));
            assert!(seq.get_signature_for_band(10).is_none());
        }

    }

//    //////////////////////////////////////////////////////
//    // MARK Sequences tests
//    //////////////////////////////////////////////////////
    mod sequences_tests {
        use super::super::*;

        #[test]
        fn test_new(){
            let in_file_name = "testInput".to_string();
            let sequences =  Sequences::new(4);
            assert_eq!(sequences.sequences.capacity(), 4);
        }
//
        #[test]
        pub fn test_add_seq(){
            let in_file_name = "testInput".to_string();
            let mut sequences =  Sequences::new(4);
            sequences.add_seq(Sequence::new(10, vec![0,1,2,3,4], 4));
            assert_eq!(sequences.sequences[0].id, 10);
        }

        #[test]
        pub fn test_load_sequences_from(){
            let (nb_sequences, kmer_size) = (3, 2);
            let file_name = "/Users/mahdi/IdeaProjects/test_fasta_readers/data/test.fa".to_string();
            let mut seqs = Sequences::new(nb_sequences);
            seqs.load_sequences_from(file_name, kmer_size, 0);
            //println!("items in sequences are {:?}", seqs.sequences);
            assert_eq!(seqs.sequences.len(), 3);
            assert_eq!(seqs.sequences[0].id, 0);
            assert_eq!(seqs.sequences[0].kmers_counts, [1, 5, 3, 10, 2, 2, 9, 5, 8, 5, 0, 17, 8, 6, 18, 0]);
            assert_eq!(seqs.sequences[1].id, 1);
            assert_eq!(seqs.sequences[1].kmers_counts, [1, 5, 3, 10, 2, 2, 9, 5, 8, 5, 0, 17, 8, 6, 18, 0]);
        }

//        #[test]
//        pub fn remove_signature_from_sequences(){
//
//            let (nb_sequences, kmer_size, sig_id, max_seqs, subset_number) = (3, 2, 1, 3, 0);
//            let file_name = "/Users/mahdi/IdeaProjects/test_fasta_readers/data/test.fa".to_string();
//            let mut seqs = Sequences::new(nb_sequences);
//            seqs.load_sequences_from(file_name, kmer_size);
//            let mut signature = Signature::new(sig_id.to_string(), max_seqs, subset_number);
//
//            signature.add_sequence(0);
//            {
//                let  seq = &mut seqs.sequences[0];
//                seq.add_signature_for_band(signature.signature_id.clone(), signature.subset_number);
//                assert_eq!(seq.signatures[0], "1".to_string());
//                assert_eq!(seq.id, signature.sequences[0]);
//            }
//
//            signature.add_sequence(2);
//            {
//                let seq = &mut seqs.sequences[2];
//                seq.add_signature_for_band(signature.signature_id.clone(), signature.subset_number);
//                assert_eq!(seq.signatures[0], "1".to_string());
//                assert_eq!(seq.id, signature.sequences[1]);
//            }
//
//            seqs.remove_signature_from_sequences(signature);
//            let ref seq = seqs.sequences[0];
//            assert_eq! (seq.signatures[0], "-1".to_string());
//            let ref seq = seqs.sequences[2];
//            assert_eq! (seq.signatures[0], "-1".to_string());
//
//        }
    }

}
