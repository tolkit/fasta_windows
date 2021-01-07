// calculate gc content/skew here

pub mod gcu8 {

    use std::collections::HashMap;

    pub struct GCStats {
        pub gc_content: f32,
        pub gc_skew: f32,
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
    
    // compute the gc content of a given sequence window.
    
    pub fn gc_content(dna: &[u8]) -> GCStats {
        // sequence length
        let length: f32 = dna.len() as f32;
        // G/C/A/T counts
        let counts = nucleotide_counts(dna);
        // G/C's
        let g_counts = counts.get(&71).unwrap_or(&0); // 71 == G
        let c_counts = counts.get(&67).unwrap_or(&0); // 67 == C
    
        GCStats {
            gc_content: ((g_counts + c_counts) as f32 / length) * 100.0,
            gc_skew: (g_counts - c_counts) as f32 / (g_counts + c_counts) as f32,
        }
    }
}