pub mod kmeru8 {

    use std::collections::HashMap;

    // still want to implement some kind of distance metric from
    // a window to the genome wide kmer profile.
    // currently canonical: true is quite a costly computation...
    pub struct KmerStats {
        pub kmer_hash: HashMap<Vec<u8>, i32>,
        pub kmer_diversity: usize,
    }

    // this function takes a chromosome (scaff/contig) of &[u8] dna
    // splits into overlapping windows of length k
    // and collects them into a hashmap
    // the length of the keys of this hashmap is the kmer diversity.

    // if canonical: true
    // the kmer and its reverse complement are considered
    // and only the lexicographically 'smaller' kmer is stored in the map.

    pub fn kmer_diversity(dna: &[u8], kmer_size: usize, canonical: bool) -> KmerStats {
        // generate all the kmers in a window
        let kmers = dna.windows(kmer_size);
        // store the kmers
        let mut map = HashMap::new();
        for kmer in kmers {
            // kmer to upper
            // unfortunately this creates a copy
            // but in place manipulation seems difficult, because rust.
            let mut kmer_upper = kmer.to_ascii_uppercase();
            if canonical {
                // switch to lexicographically lower kmer
                let rev_kmer = reverse_complement(&kmer_upper);
                if rev_kmer < kmer_upper {
                    kmer_upper = rev_kmer;
                }
                // skip where kmer contains an N (or any other invalid character?)
                if kmer_upper.contains(&b'N') {
                    continue;
                }
                let count = map.entry(kmer_upper).or_insert(0);
                *count += 1;
            } else {
                if kmer_upper.contains(&b'N') {
                    continue;
                }
                let count = map.entry(kmer_upper).or_insert(0);
                *count += 1;
            }
        }
        let diversity = map.keys().len();

        KmerStats {
            kmer_hash: map,
            // do something with the hashmap
            kmer_diversity: diversity,
        }
    }

    fn reverse_complement(dna: &[u8]) -> Vec<u8> {
        let dna_vec = dna.to_vec();
        let mut revcomp = Vec::new();

        for base in dna_vec.iter() {
            revcomp.push(switch_base(*base))
        }
        revcomp.as_mut_slice().reverse();
        revcomp
    }

    // works on uppercase ascii, so
    // no need for lowercase here.

    fn switch_base(c: u8) -> u8 {
        match c {
            b'A' => b'T',
            b'C' => b'G',
            b'T' => b'A',
            b'G' => b'C',
            b'N' => b'N',
            _ => b'N',
        }
    }
}
