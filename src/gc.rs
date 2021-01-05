// calculate gc content/skew here

pub mod gc {

use std::collections::HashMap;

// function below reveals other ambiguous bases present in assemblies, not sure
// how to deal with those yet.

// count the number of each nucleotide in a given sequence window.

pub fn nucleotide_counts(dna: &str) -> HashMap<char, i32> {
    let mut map = HashMap::new();

    for nucleotide in dna.chars() {
        let count = map.entry(nucleotide).or_insert(0);
        *count += 1;
    }
    map
}

pub struct GCStats {
    pub gc_content: f32,
    pub gc_skew: f32,
}

// compute the gc content of a given sequence window.

pub fn gc_content(dna: &str) -> GCStats {
    // sequence length
    let length: f32 = dna.len() as f32;
    // G/C/A/T counts
    let counts = nucleotide_counts(dna);
    // G/C's
    let g_counts = counts.get(&'G').unwrap_or(&0);
    let c_counts = counts.get(&'C').unwrap_or(&0);

    GCStats {
        gc_content: ((g_counts + c_counts) as f32 / length) * 100.0,
        gc_skew: (g_counts - c_counts) as f32 / (g_counts + c_counts) as f32,
    }
}

}