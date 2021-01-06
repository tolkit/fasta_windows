pub mod kmer {

    use std::collections::HashMap;

    use crate::windows::windows;
    
    // still want to implement some kind of distance metric from 
    // a window to the genome wide kmer profile.

    // currently canonical: true is quite a costly computation...

    pub struct KmerStats {
        pub kmer_hash: HashMap<String, i32>,
        pub kmer_diversity: usize,
    }

    pub fn kmer_diversity(dna: &str, kmer_size: usize, canonical: bool) -> KmerStats {
        // generate all the kmers in a window
        let kmers = windows::char_windows(&dna, kmer_size, 1);
        
        let mut map = HashMap::new();
        
        for mut kmer in kmers {
            if canonical {
                // switch to lexicographically lower kmer
                let rev_kmer = revcomp(&kmer);
                if rev_kmer < kmer.to_string() {
                    kmer = &rev_kmer;
                }
                // skip where kmer contains an N
                if kmer.contains("N"){ 
                    continue; 
                } 
                let count = map.entry(kmer.to_string()).or_insert(0);
                *count += 1;
            } else {
                if kmer.contains("N"){ 
                    continue; 
                }
                let count = map.entry(kmer.to_string()).or_insert(0);
                *count += 1;
            }
        }
        let diversity = map.keys().len();

        KmerStats {
            kmer_hash: map,
            kmer_diversity: diversity,
        }
    }

// this is sloooow.
// reverse-complement of a dna sequence
fn revcomp(dna: &str) -> String {
        let mut rc_dna: String = String::with_capacity(dna.len());
        for c in dna.chars().rev() {
            rc_dna.push(switch_base(c))
        }
        rc_dna
    }

// Complementary nucleotide
fn switch_base(c: char) -> char {
        match c {
            'A' => 'T',
            'C' => 'G',
            'T' => 'A',
            'G' => 'C',
            'N' => 'N',
            // anything else, N
            _ => 'N', 
        }
    }
}