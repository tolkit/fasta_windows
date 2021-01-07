pub mod kmeru8 {

    use std::collections::HashMap;
    
    // still want to implement some kind of distance metric from 
    // a window to the genome wide kmer profile.
    
    // currently canonical: true is quite a costly computation...
    pub struct KmerStats {
        pub kmer_hash: HashMap<Vec<u8>, i32>,
        pub kmer_diversity: usize,
    }

    pub fn kmer_diversity(dna: &[u8], kmer_size: usize, canonical: bool) -> KmerStats {
        // generate all the kmers in a window
        let kmers = dna.windows(kmer_size);
        // store the kmers
        let mut map = HashMap::new();
        
        for mut kmer in kmers {
            if canonical {
                // switch to lexicographically lower kmer
                // probably inefficient...
                let rev_kmer = reverse_complement(kmer);
                if rev_kmer < kmer.to_vec() {
                    kmer = &rev_kmer;
                }
                // skip where kmer contains an N
                if kmer.contains(&b'N'){ 
                    continue; 
                } 
                let count = map.entry(kmer.to_vec()).or_insert(0);
                *count += 1;
            } else {
                if kmer.contains(&b'N'){ 
                    continue; 
                }
                let count = map.entry(kmer.to_vec()).or_insert(0);
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