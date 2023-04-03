use crate::kmer_maps::KmerMap;
use rayon::prelude::*;
use std::sync::mpsc::channel;

// calculating shannon diversity of di/tri/tetranucleotides
// convenience struct
pub struct KmerStats {
    // the frequency distribution of the kmer
    pub freq_dist_k: Vec<i32>,
    // the length of the kmer
    pub kmer_length: usize,
    // the shannon diversity
    pub shannon: f64,
}

// punted to main.rs
pub struct ShannonDiversity {
    // of dinucleotides
    pub dinucleotides: f64,
    // of trinucleotides
    pub trinucleotides: f64,
    // of tetranucleotides
    pub tetranucleotides: f64,
    // the dimer frequency distribution
    pub di_freq: Vec<i32>,
    // the trimer frequency distribution
    pub tri_freq: Vec<i32>,
    // the tetramer frequency distribution
    pub tetra_freq: Vec<i32>,
}

pub fn kmer_diversity(dna: &[u8], kmer_maps: Vec<KmerMap>) -> ShannonDiversity {
    // parallel iterate over 2-4mers
    let (sender, receiver) = channel();
    kmer_maps.into_par_iter().for_each_with(sender, |s, i| {
        // need to do di/tri/tetranucleotides
        // generate all the kmers in a window
        let kmers = dna.windows(i.len);
        let mut map = i.map;

        // iterate over sliding windows of length k
        for kmer in kmers {
            // kmer to upper
            // unfortunately this creates a copy
            // but in place manipulation seems difficult, because rust.
            let kmer_upper = kmer.to_ascii_uppercase();
            if kmer_upper.contains(&b'N') {
                continue;
            }
            let count = map.entry(kmer_upper).or_insert(0);
            *count += 1;
        }
        // now calculate shannon diversity
        // is this right??
        let shannon = shannon_diversity(map.values().cloned().collect());

        // save the hashmap here too.
        // HashMap -> Vec -> sorted Vec by HashMap keys
        // so keys should always be in the same order.
        let mut map_vec: Vec<_> = map.into_iter().collect();
        map_vec.sort_by(|x, y| x.0.cmp(&y.0));
        let values: Vec<i32> = map_vec.iter().map(|(_x, y)| *y).collect();

        s.send(KmerStats {
            freq_dist_k: values,
            kmer_length: i.len,
            shannon,
        })
        .expect("KmerStats did not send!");
    });
    // collect stats
    let kmer_stats: Vec<KmerStats> = receiver.iter().collect();

    // decompose into separate shannon indices
    // and the k-mer freq spectra
    // TODO: is there a better way to do this?
    let mut dinucleotides: f64 = 0.0;
    let mut trinucleotides: f64 = 0.0;
    let mut tetranucleotides: f64 = 0.0;
    let mut divalues: Vec<i32> = Vec::new();
    let mut trivalues: Vec<i32> = Vec::new();
    let mut tetravalues: Vec<i32> = Vec::new();

    for stat in kmer_stats {
        match stat.kmer_length {
            2usize => {
                dinucleotides += stat.shannon;
                divalues = stat.freq_dist_k;
            }
            3usize => {
                trinucleotides += stat.shannon;
                trivalues = stat.freq_dist_k;
            }
            4usize => {
                tetranucleotides += stat.shannon;
                tetravalues = stat.freq_dist_k;
            }
            _ => (),
        }
    }

    ShannonDiversity {
        dinucleotides,
        trinucleotides,
        tetranucleotides,
        di_freq: divalues,
        tri_freq: trivalues,
        tetra_freq: tetravalues,
    }
}

// using the natural log
fn shannon_diversity(vec: Vec<i32>) -> f64 {
    // sum elements to get proportions
    let vec_sum: i32 = vec.iter().sum();
    let mut diversity = 0f64;

    for count in vec.iter().filter(|count| **count > 0i32) {
        let probability = *count as f64 / (vec_sum as f64);
        diversity -= probability * probability.log2();
    }
    diversity
}

pub fn reverse_complement(dna: &[u8]) -> Vec<u8> {
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

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_reverse_complement() {
        let short_dna_string = "AACCTTGG".as_bytes();

        let revcomp_dna = reverse_complement(short_dna_string);

        let should_be = "CCAAGGTT".as_bytes().to_vec();

        assert_eq!(should_be, revcomp_dna);
    }

    #[test]
    fn test_shannon_entropy() {
        todo!()
    }
}
