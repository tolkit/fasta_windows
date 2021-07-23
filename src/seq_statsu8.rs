// calculate gc content/skew here

pub mod seq_statsu8 {

    use std::collections::HashMap;

    // TODO: add AT skew? AT content?
    pub struct SeqStats {
        pub gc_proportion: f32,
        pub gc_skew: f32,
        pub shannon_entropy: f64,
        // proportions of nucleotides (& N's)
        pub g_s: f32,
        pub c_s: f32,
        pub a_s: f32,
        pub t_s: f32,
        pub n_s: f32,
    }
    // function below reveals other ambiguous bases present in assemblies, not sure
    // how to deal with those yet.
    // count the number of each nucleotide in a given sequence window.
    // I've learned that bio::fasta does UTF8 checks, so we can be confident
    // of all the characters in the sequence
    fn nucleotide_counts(dna: &[u8]) -> HashMap<&u8, i32> {
        let mut map = HashMap::new();
        for nucleotide in dna {
            let count = map.entry(nucleotide).or_insert(0);
            *count += 1;
        }
        map
    }

    // store each u8 in a hashmap, extract the values and summarise
    // add in here option to pass over lowercase letters
    pub fn seq_stats(dna: &[u8], masked: bool) -> SeqStats {
        // sequence length
        let length: f32 = dna.len() as f32;
        // G/C/A/T counts
        let counts = nucleotide_counts(dna);
        // upper and lower cases accounted for.

        let g_counts: i32;
        let c_counts: i32;
        let a_counts: i32;
        let t_counts: i32;
        let n_counts: i32;

        if masked {
            g_counts = *counts.get(&71).unwrap_or(&0); // 71 == G;
            c_counts = *counts.get(&67).unwrap_or(&0); // 67 == C;
            a_counts = *counts.get(&65).unwrap_or(&0); // 65 == A;
            t_counts = *counts.get(&84).unwrap_or(&0); // 84 == T;
            n_counts = *counts.get(&78).unwrap_or(&0); // 78 == N;
        } else {
            g_counts = counts.get(&71).unwrap_or(&0) + counts.get(&103).unwrap_or(&0); // 71 == G; 103 == g
            c_counts = counts.get(&67).unwrap_or(&0) + counts.get(&99).unwrap_or(&0); // 67 == C; 99 == c
            a_counts = counts.get(&65).unwrap_or(&0) + counts.get(&97).unwrap_or(&0); // 65 == A; 97 == a
            t_counts = counts.get(&84).unwrap_or(&0) + counts.get(&116).unwrap_or(&0); // 84 == T; 116 == t
            n_counts = counts.get(&78).unwrap_or(&0) + counts.get(&110).unwrap_or(&0);
            // 78 == N; 110 == n
        }

        // shannon entropy of the window
        // see https://github.com/fkie-cad/entropython/blob/main/src/lib.rs

        let mut byte_count = [0u64; 256];
        for byte in dna {
            // change lowercase nucleotides to uppercase ones.
            match byte {
                103u8 => byte_count[71u8 as usize] += 1,
                99u8 => byte_count[67u8 as usize] += 1,
                97u8 => byte_count[65u8 as usize] += 1,
                116u8 => byte_count[84u8 as usize] += 1,
                110u8 => byte_count[78u8 as usize] += 1,
                _ => byte_count[*byte as usize] += 1,
            }
        }
        let mut entropy = 0f64;
        for counted_num in byte_count.iter().filter(|num| **num > 0u64) {
            let byte_probability = *counted_num as f64 / (dna.len() as f64);
            entropy -= byte_probability * byte_probability.log2();
        }
        SeqStats {
            gc_proportion: ((g_counts + c_counts) as f32
                / (g_counts + c_counts + a_counts + t_counts) as f32),
            gc_skew: (g_counts - c_counts) as f32 / (g_counts + c_counts) as f32,
            shannon_entropy: entropy,
            g_s: ((g_counts) as f32 / length),
            c_s: ((c_counts) as f32 / length),
            a_s: ((a_counts) as f32 / length),
            t_s: ((t_counts) as f32 / length),
            n_s: ((n_counts) as f32 / length),
        }
    }
}
