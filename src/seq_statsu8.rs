// calculate gc content/skew here

pub mod seq_statsu8 {

    use std::collections::HashMap;

    // TODO: add at skew? at content?
    pub struct SeqStats {
        pub gc_content: f32,
        pub gc_skew: f32,
        pub shannon_entropy: f64,
        pub g_s: i32,
        pub c_s: i32,
        pub a_s: i32,
        pub t_s: i32,
        pub n_s: i32,
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

    pub fn seq_stats(dna: &[u8]) -> SeqStats {
        // sequence length
        let length: f32 = dna.len() as f32;
        // G/C/A/T counts
        let counts = nucleotide_counts(dna);
        // upper and lower cases accounted for.
        let g_counts = counts.get(&71).unwrap_or(&0) + counts.get(&103).unwrap_or(&0); // 71 == G; 103 == g
        let c_counts = counts.get(&67).unwrap_or(&0) + counts.get(&99).unwrap_or(&0); // 67 == C; 99 == c
        let a_counts = counts.get(&65).unwrap_or(&0) + counts.get(&97).unwrap_or(&0); // 65 == A; 97 == a
        let t_counts = counts.get(&84).unwrap_or(&0) + counts.get(&116).unwrap_or(&0); // 84 == T; 116 == t
        let n_counts = counts.get(&78).unwrap_or(&0) + counts.get(&110).unwrap_or(&0); // 78 == N; 110 == n

        // shannon entropy of the window
        // see https://github.com/fkie-cad/entropython/blob/main/src/lib.rs
        // TODO: will this be higher if lower/uppercases are not synonymised?
        let mut byte_count = [0u64; 256];
        for byte in dna {
            byte_count[*byte as usize] += 1;
        }
        let mut entropy = 0f64;
        for counted_num in byte_count.iter().filter(|num| **num > 0u64) {
            let byte_probability = *counted_num as f64 / (dna.len() as f64);
            entropy -= byte_probability * byte_probability.log2();
        }
        SeqStats {
            gc_content: ((g_counts + c_counts) as f32 / length) * 100.0,
            gc_skew: (g_counts - c_counts) as f32 / (g_counts + c_counts) as f32,
            shannon_entropy: entropy,
            g_s: g_counts,
            c_s: c_counts,
            a_s: a_counts,
            t_s: t_counts,
            n_s: n_counts,
        }
    }
}
