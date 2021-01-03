pub mod kmer {
    use std::collections::HashMap;

    use crate::windows::windows;

    pub fn kmer_diversity(dna: &str, kmer_size: usize) -> usize {
        // generate all the kmers in a window
        let kmers = windows::char_windows(dna, kmer_size, 1);

        let mut map = HashMap::new();
        for kmer in kmers {
                let count = map.entry(kmer).or_insert(0);
                *count += 1;
            }
        // simply return the length of the keys
        map.keys().len()
        // TODO compare the kmer frequency spectrum with the genome wide signature.
        // see dnaglider on tips
    }

}