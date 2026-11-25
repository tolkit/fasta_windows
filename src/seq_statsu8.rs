pub struct SeqStats {
    pub gc_proportion: f32,
    pub gc_skew: f32,
    pub at_skew: f32,
    pub shannon_entropy: f64,
    /// counts of A, C, G, T, N
    pub nuc_counts: [u64; 5],
    // proportions of nucleotides (& N's)
    pub g_s: f32,
    pub c_s: f32,
    pub a_s: f32,
    pub t_s: f32,
    pub n_s: f32,
    pub masked: f32,
    pub len: f32,
}

/// Count every byte value in `dna`.
///
/// Four interleaved tables break the store-to-load dependency that a single
/// table hits on runs of the same base (homopolymers, satellites). The u32
/// tables are flushed into u64 totals every 2^30 bytes so very large windows
/// cannot overflow.
pub fn byte_histogram(dna: &[u8]) -> [u64; 256] {
    let mut total = [0u64; 256];
    for block in dna.chunks(1 << 30) {
        let mut h = [[0u32; 256]; 4];
        let mut quads = block.chunks_exact(4);
        for q in &mut quads {
            h[0][q[0] as usize] += 1;
            h[1][q[1] as usize] += 1;
            h[2][q[2] as usize] += 1;
            h[3][q[3] as usize] += 1;
        }
        for &b in quads.remainder() {
            h[0][b as usize] += 1;
        }
        for i in 0..256 {
            total[i] += (h[0][i] + h[1][i] + h[2][i] + h[3][i]) as u64;
        }
    }
    total
}

/// All base-composition statistics of a window, from a single pass over it.
pub fn seq_stats(dna: &[u8], masked: bool) -> SeqStats {
    let length = dna.len() as f32;
    let hist = byte_histogram(dna);
    let n = |b: u8| hist[b as usize] as i64;

    let (g, c, a, t, n_count, masked_count, w, s);
    if masked {
        // only upper case bases count
        g = n(b'G');
        c = n(b'C');
        a = n(b'A');
        t = n(b'T');
        n_count = n(b'N');
        masked_count = 0;
        w = n(b'W');
        s = n(b'S');
    } else {
        g = n(b'G') + n(b'g');
        c = n(b'C') + n(b'c');
        a = n(b'A') + n(b'a');
        t = n(b'T') + n(b't');
        n_count = n(b'N') + n(b'n');
        // all valid lower case bases
        masked_count = b"acgtmrwsykvhbdn".iter().map(|&b| n(b)).sum();
        // additional weak (ambiguous A/T) and strong (ambiguous G/C) bases
        w = n(b'W') + n(b'w');
        s = n(b'S') + n(b's');
    }

    // Shannon entropy over byte values, with lower case acgtn folded into
    // upper case; other bytes (including lower case IUPAC codes) keep their
    // own bin. See https://github.com/fkie-cad/entropython/blob/main/src/lib.rs
    let mut bins = hist;
    for (lower, upper) in [
        (b'g', b'G'),
        (b'c', b'C'),
        (b'a', b'A'),
        (b't', b'T'),
        (b'n', b'N'),
    ] {
        bins[upper as usize] += bins[lower as usize];
        bins[lower as usize] = 0;
    }
    let mut entropy = 0f64;
    for &count in bins.iter().filter(|&&count| count > 0) {
        let p = count as f64 / dna.len() as f64;
        entropy -= p * p.log2();
    }

    SeqStats {
        gc_proportion: (g + c + s) as f32 / (g + c + s + a + t + w) as f32,
        gc_skew: (g - c) as f32 / (g + c) as f32,
        at_skew: (a - t) as f32 / (a + t) as f32,
        shannon_entropy: entropy,
        nuc_counts: [a as u64, c as u64, g as u64, t as u64, n_count as u64],
        g_s: g as f32 / length,
        c_s: c as f32 / length,
        a_s: a as f32 / length,
        t_s: t as f32 / length,
        n_s: n_count as f32 / length,
        masked: masked_count as f32 / length,
        len: length,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_byte_histogram() {
        let dna = b"AACCTTGGa";
        let h = byte_histogram(dna);
        assert_eq!(h[b'A' as usize], 2);
        assert_eq!(h[b'C' as usize], 2);
        assert_eq!(h[b'G' as usize], 2);
        assert_eq!(h[b'T' as usize], 2);
        assert_eq!(h[b'a' as usize], 1);
        assert_eq!(h.iter().sum::<u64>(), dna.len() as u64);
        // lengths that are not a multiple of the 4-way unroll
        for len in 0..9 {
            assert_eq!(byte_histogram(&dna[..len]).iter().sum::<u64>(), len as u64);
        }
    }

    #[test]
    fn test_nucleotide_counts() {
        let stats = seq_stats(b"AACCTTGGNacgtn", false);
        assert_eq!(stats.nuc_counts, [3, 3, 3, 3, 2]);
        let stats = seq_stats(b"AACCTTGGNacgtn", true);
        assert_eq!(stats.nuc_counts, [2, 2, 2, 2, 1]);
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

    #[test]
    fn test_entropy_folds_case() {
        let upper = seq_stats(b"ACGTACGT", false).shannon_entropy;
        let mixed = seq_stats(b"acgtACGT", false).shannon_entropy;
        assert_eq!(upper, 2.0);
        assert_eq!(upper, mixed);
    }
}
