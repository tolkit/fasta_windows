// calculate gc content/skew here

pub mod seq_statsu8 {

    use std::collections::HashMap;

    pub struct SeqStats {
        pub gc_content: f32,
        pub gc_skew: f32,
        pub shannon_entropy: f64,
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
    
    pub fn seq_stats(dna: &[u8]) -> SeqStats {
        // sequence length
        let length: f32 = dna.len() as f32;
        // G/C/A/T counts
        let counts = nucleotide_counts(dna);
        // G/C's
        let g_counts = counts.get(&71).unwrap_or(&0); // 71 == G
        let c_counts = counts.get(&67).unwrap_or(&0); // 67 == C

        // shannon entropy of the window
        // see https://github.com/fkie-cad/entropython/blob/main/src/lib.rs
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
        }
    }
}