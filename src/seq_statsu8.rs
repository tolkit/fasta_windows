use std::collections::HashMap;

pub struct SeqStats {
    pub gc_proportion: f32,
    pub gc_skew: f32,
    pub at_skew: f32,
    pub shannon_entropy: f64,
    pub nuc_counts: Vec<i32>,
    // proportions of nucleotides (& N's)
    pub g_s: f32,
    pub c_s: f32,
    pub a_s: f32,
    pub t_s: f32,
    pub n_s: f32,
    pub masked: f32,
    pub len: f32,
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
    let masked_counts: i32;
    let w_counts: i32;
    let s_counts: i32;

    if masked {
        g_counts = *counts.get(&b'G').unwrap_or(&0);
        c_counts = *counts.get(&b'C').unwrap_or(&0);
        a_counts = *counts.get(&b'A').unwrap_or(&0);
        t_counts = *counts.get(&b'T').unwrap_or(&0);
        n_counts = *counts.get(&b'N').unwrap_or(&0);
        masked_counts = 0;
        w_counts = *counts.get(&b'W').unwrap_or(&0);
        s_counts = *counts.get(&b'S').unwrap_or(&0);
    } else {
        g_counts = counts.get(&b'G').unwrap_or(&0) + counts.get(&b'g').unwrap_or(&0);
        c_counts = counts.get(&b'C').unwrap_or(&0) + counts.get(&b'c').unwrap_or(&0);
        a_counts = counts.get(&b'A').unwrap_or(&0) + counts.get(&b'a').unwrap_or(&0);
        t_counts = counts.get(&b'T').unwrap_or(&0) + counts.get(&b't').unwrap_or(&0);
        n_counts = counts.get(&b'N').unwrap_or(&0) + counts.get(&b'n').unwrap_or(&0);
        // All valid lower case bases
        masked_counts = counts.get(&b'a').unwrap_or(&0)
            + counts.get(&b'c').unwrap_or(&0)
            + counts.get(&b'g').unwrap_or(&0)
            + counts.get(&b't').unwrap_or(&0)
            + counts.get(&b'm').unwrap_or(&0)
            + counts.get(&b'r').unwrap_or(&0)
            + counts.get(&b'w').unwrap_or(&0)
            + counts.get(&b's').unwrap_or(&0)
            + counts.get(&b'y').unwrap_or(&0)
            + counts.get(&b'k').unwrap_or(&0)
            + counts.get(&b'v').unwrap_or(&0)
            + counts.get(&b'h').unwrap_or(&0)
            + counts.get(&b'b').unwrap_or(&0)
            + counts.get(&b'd').unwrap_or(&0)
            + counts.get(&b'n').unwrap_or(&0);
        // additional weak bases (ambiguous A/T)
        w_counts = counts.get(&b'W').unwrap_or(&0) + counts.get(&b'w').unwrap_or(&0);
        // additional strong bases (ambiguous G/C)
        s_counts = counts.get(&b'S').unwrap_or(&0) + counts.get(&b's').unwrap_or(&0);
    }

    // shannon entropy of the window
    // see https://github.com/fkie-cad/entropython/blob/main/src/lib.rs

    let mut byte_count = [0u64; 256];
    for byte in dna {
        // change lowercase nucleotides to uppercase ones.
        match byte {
            b'g' => byte_count[b'G' as usize] += 1,
            b'c' => byte_count[b'C' as usize] += 1,
            b'a' => byte_count[b'A' as usize] += 1,
            b't' => byte_count[b'T' as usize] += 1,
            b'n' => byte_count[b'N' as usize] += 1,
            _ => byte_count[*byte as usize] += 1,
        }
    }
    let mut entropy = 0f64;
    for counted_num in byte_count.iter().filter(|num| **num > 0u64) {
        let byte_probability = *counted_num as f64 / (dna.len() as f64);
        entropy -= byte_probability * byte_probability.log2();
    }
    SeqStats {
        gc_proportion: ((g_counts + c_counts + s_counts) as f32
            / (g_counts + c_counts + s_counts + a_counts + t_counts + w_counts) as f32),
        gc_skew: (g_counts - c_counts) as f32 / (g_counts + c_counts) as f32,
        at_skew: (a_counts - t_counts) as f32 / (a_counts + t_counts) as f32,
        shannon_entropy: entropy,
        nuc_counts: vec![a_counts, c_counts, g_counts, t_counts, n_counts],
        g_s: ((g_counts) as f32 / length),
        c_s: ((c_counts) as f32 / length),
        a_s: ((a_counts) as f32 / length),
        t_s: ((t_counts) as f32 / length),
        n_s: ((n_counts) as f32 / length),
        masked: ((masked_counts) as f32 / length),
        len: length,
    }
}

#[cfg(test)]
mod tests {

    use crate::seq_statsu8::seq_stats;

    use super::nucleotide_counts;

    const A: u8 = b'A';
    const C: u8 = b'C';
    const G: u8 = b'G';
    const T: u8 = b'T';

    #[test]
    fn test_nucleotide_counts() {
        let short_dna_string = "AACCTTGG".as_bytes();

        let nuc_counts = nucleotide_counts(short_dna_string);

        // Two of each!
        assert_eq!(2, *nuc_counts.get(&A).unwrap());
        assert_eq!(2, *nuc_counts.get(&C).unwrap());
        assert_eq!(2, *nuc_counts.get(&G).unwrap());
        assert_eq!(2, *nuc_counts.get(&T).unwrap());
    }

    #[test]
    fn test_masked_proportion() {
        let short_dna_string = "AAaCCcTTtGGg".as_bytes();

        let stats = seq_stats(short_dna_string, false);

        // One third masked
        assert_eq!(1.0 / 3.0, stats.masked);
    }

    #[test]
    fn test_ambiguous_gc_proportion() {
        let short_dna_string = "AASCTTGsWw".as_bytes();

        let stats = seq_stats(short_dna_string, false);

        // 4 out of 10 bases are GC
        assert_eq!(0.4, stats.gc_proportion);
    }
}
