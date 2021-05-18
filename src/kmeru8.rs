pub mod kmeru8 {

    use rayon::prelude::*;
    use std::collections::HashMap;
    use std::sync::mpsc::channel;

    // calculating shannon diversity of di/tri/tetranucleotides

    // convenience struct
    pub struct KmerStats {
        pub kmer_length: i32,
        pub shannon: f64,
    }

    // punted to main.rs
    pub struct ShannonDiversity {
        pub dinucleotides: f64,
        pub trinucleotides: f64,
        pub tetranucleotides: f64,
    }

    // if canonical: true
    // the kmer and its reverse complement are considered
    // and only the lexicographically 'smaller' kmer is stored in the map.
    // I think this should still work given the restructure? check...

    pub fn kmer_diversity(dna: &[u8], canonical: bool) -> ShannonDiversity {
        // parallel iterate over 2-4mers
        let (sender, receiver) = channel();
        (2usize..5usize)
            .into_par_iter()
            .for_each_with(sender, |s, i| {
                // need to do di/tri/tetranucleotides
                // generate all the kmers in a window
                let kmers = dna.windows(i);
                // store the kmers
                let kmer_map = gen_all_kmers(i);
                let kmers_str: Vec<&str> = kmer_map.iter().map(|s| &**s).collect();
                let mut kmers_u8: Vec<Vec<u8>> = Vec::new();

                for i in kmers_str {
                    kmers_u8.push(i.as_bytes().to_vec()); // to vec copies...
                }

                // initiate a map
                let mut map = HashMap::new();
                // insert all the possible kmers as keys
                for k in kmers_u8 {
                    map.insert(k, 0i32);
                }

                // iterate over sliding windows of length k
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
                        if kmer_upper.contains(&b'N') || kmer_upper.contains(&b'n') {
                            continue;
                        }
                        let count = map.entry(kmer_upper).or_insert(0);
                        *count += 1;
                    } else {
                        if kmer_upper.contains(&b'N') || kmer_upper.contains(&b'n') {
                            continue;
                        }
                        let count = map.entry(kmer_upper).or_insert(0);
                        *count += 1;
                    }
                }
                // now calculate shannon diversity
                let shannon = shannon_diversity(map.values().cloned().collect());

                s.send(KmerStats {
                    kmer_length: i as i32,
                    shannon: shannon,
                })
                .expect("KmerStats did not send!");
            });
        // collect stats
        let kmer_stats: Vec<KmerStats> = receiver.iter().collect();

        // decompose into separate shannon indices
        let mut dinucleotides: f64 = 0.0;
        let mut trinucleotides: f64 = 0.0;
        let mut tetranucleotides: f64 = 0.0;

        for stat in kmer_stats {
            match stat.kmer_length {
                2i32 => dinucleotides += stat.shannon,
                3i32 => trinucleotides += stat.shannon,
                4i32 => tetranucleotides += stat.shannon,
                _ => (),
            }
        }

        ShannonDiversity {
            dinucleotides: dinucleotides,
            trinucleotides: trinucleotides,
            tetranucleotides: tetranucleotides,
        }
    }

    // Nice recursive function from
    // https://github.com/szhongren/b363-bioinformatics/blob/2495b36575887222f5b1c266046fdd6ded23bf68/lab/challenge/src/main.rs
    fn gen_all_kmers(k: usize) -> Vec<String> {
        let mut output: Vec<String> = Vec::new();
        let nucleotides = ["A", "C", "G", "T"].iter().map(|x| x.to_string()).collect();
        if k == 0 {
            return output;
        } else if k == 1 {
            return nucleotides;
        } else {
            for suffix in gen_all_kmers(k - 1) {
                for nt in &nucleotides {
                    let new = suffix.clone() + &nt;
                    output.push(new);
                }
            }
        }
        output
    }

    // using the natural log
    fn shannon_diversity(vec: Vec<i32>) -> f64 {
        // sum elements to get proportions
        let vec_sum: i32 = vec.iter().sum();
        // collect the proportions * ln(proportions)
        let mut prop_logs: Vec<f64> = Vec::new();

        for e in vec {
            if e > 0 {
                prop_logs.push((e as f64 / vec_sum as f64) * ((e as f64) / vec_sum as f64).ln());
            }
        }

        // take the sum of these proportions * ln(proportions)
        let res: f64 = prop_logs.iter().sum();
        res * -1.0
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
