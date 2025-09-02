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

/// Zero-order KT code length for DNA (A,C,G,T), returned as bits per base.
/// Skips non-ACGT.
fn kt0_bits_per_base_dna(seq: &[u8]) -> f64 {
    #[inline]
    fn nuc_to_sym(b: u8) -> Option<usize> {
        match b {
            b'A' | b'a' => Some(0),
            b'C' | b'c' => Some(1),
            b'G' | b'g' => Some(2),
            b'T' | b't' => Some(3),
            _ => None,
        }
    }
    let m = 4usize;
    let mut counts = [0u32; 4];
    let mut n_eff = 0usize;
    let mut sum_log2 = 0.0;

    for &b in seq {
        if let Some(sym) = nuc_to_sym(b) {
            let c_s = counts[sym] as f64;
            let n = (counts[0] + counts[1] + counts[2] + counts[3]) as f64;
            let num = c_s + 0.5;
            let den = n + (m as f64) / 2.0; // N + m/2  (m=4 => N+2)
            sum_log2 += (num / den).ln() / std::f64::consts::LN_2;
            counts[sym] += 1;
            n_eff += 1;
        }
    }
    if n_eff == 0 {
        return 0.0;
    }
    let bits = -sum_log2;
    bits / (n_eff as f64)
}

// Context Tree Weighting (CTW) for DNA (A,C,G,T).
// Returns bits per effective base for the given window using max context depth `max_depth`.
// Non-ACGT symbols are skipped and flush the context.
//
// References:
// - Willems, Shtarkov & Tjalkens (1995): CTW basics.
// - Multinomial KT estimator with 1/2 pseudo-counts.
// - m-ary CTW mixture weight beta = 1 - 2^{-(m-1)} (=> 7/8 for DNA).
/// NOTE: If `max_depth == 0` this returns the exact KT(0) code length.
pub fn ctw_bits_per_base_dna(dna: &[u8], max_depth: usize) -> f64 {
    if max_depth == 0 {
        return kt0_bits_per_base_dna(dna);
    }

    const M: usize = 4; // A,C,G,T
    let beta: f64 = 0.5; // standard CTW mixture

    #[inline]
    fn nuc_to_sym(b: u8) -> Option<usize> {
        match b {
            b'A' | b'a' => Some(0),
            b'C' | b'c' => Some(1),
            b'G' | b'g' => Some(2),
            b'T' | b't' => Some(3),
            _ => None,
        }
    }

    #[inline]
    fn log2(x: f64) -> f64 {
        x.ln() / std::f64::consts::LN_2
    }

    #[inline]
    fn log2_sum_weighted(a_log2: f64, b_log2: f64, beta: f64) -> f64 {
        // log2( beta*2^a + (1-beta)*2^b )
        if !a_log2.is_finite() && !b_log2.is_finite() {
            return f64::NEG_INFINITY;
        }
        let m = a_log2.max(b_log2);
        let ta = if (a_log2 - m) < -50.0 {
            0.0
        } else {
            beta * (2f64).powf(a_log2 - m)
        };
        let tb = if (b_log2 - m) < -50.0 {
            0.0
        } else {
            (1.0 - beta) * (2f64).powf(b_log2 - m)
        };
        m + log2(ta + tb)
    }

    struct Node {
        m: usize,
        beta: f64,
        counts: Vec<u32>,
        total: u32,
        log_p_kt: f64,                    // log2 P_KT at this node
        log_w: f64,                       // log2 weighted prob at this node
        children: Vec<Option<Box<Node>>>, // arity m (4)
    }

    impl Node {
        fn new(m: usize, beta: f64) -> Self {
            Self {
                m,
                beta,
                counts: vec![0; m],
                total: 0,
                log_p_kt: 0.0, // prob 1
                log_w: 0.0,    // prob 1
                children: (0..m).map(|_| None).collect(),
            }
        }

        #[inline]
        fn ensure_child(&mut self, a: usize) -> &mut Node {
            if self.children[a].is_none() {
                self.children[a] = Some(Box::new(Node::new(self.m, self.beta)));
            }
            self.children[a].as_mut().unwrap()
        }

        /// Update with `sym` given `ctx` (most-recent-first).
        /// Leaf rule: if ctx.is_empty(), log_w := log_p_kt (no mixture).
        fn update(&mut self, ctx: &[usize], sym: usize) {
            if let Some((&a, rest)) = ctx.split_first() {
                // 1) deeper
                self.ensure_child(a).update(rest, sym);

                // 2) KT update
                let c_s_old = self.counts[sym] as f64;
                let n_old = self.total as f64;
                let num = c_s_old + 0.5;
                let den = n_old + (self.m as f64) / 2.0;
                self.log_p_kt += (num / den).ln() / std::f64::consts::LN_2;

                self.counts[sym] += 1;
                self.total += 1;

                // 3) children product (sum of logs)
                let mut sum_children_log_w = 0.0;
                for ch in self.children.iter() {
                    if let Some(ref child) = ch {
                        sum_children_log_w += child.log_w;
                    }
                }

                // 4) mixture
                self.log_w = log2_sum_weighted(self.log_p_kt, sum_children_log_w, self.beta);
            } else {
                // leaf
                let c_s_old = self.counts[sym] as f64;
                let n_old = self.total as f64;
                let num = c_s_old + 0.5;
                let den = n_old + (self.m as f64) / 2.0;
                self.log_p_kt += (num / den).ln() / std::f64::consts::LN_2;

                self.counts[sym] += 1;
                self.total += 1;

                self.log_w = self.log_p_kt; // no mixture at leaves
            }
        }
    }

    // Online: sum_t log2 (W_after / W_before) at root
    let mut root = Node::new(M, beta);
    let mut ctx_vec: Vec<usize> = Vec::with_capacity(max_depth);

    let mut total_delta_logw = 0.0;
    let mut n_eff = 0usize;

    for &b in dna {
        let Some(sym) = nuc_to_sym(b) else {
            ctx_vec.clear(); // flush on N, masked, etc.
            continue;
        };

        let before = root.log_w;
        root.update(&ctx_vec, sym);
        let after = root.log_w;

        total_delta_logw += after - before;
        n_eff += 1;

        if ctx_vec.len() == max_depth {
            ctx_vec.pop();
        }
        ctx_vec.insert(0, sym);
    }

    if n_eff == 0 {
        return 0.0;
    }
    let bits = -total_delta_logw;
    bits / (n_eff as f64)
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

    // --- CTW tests ---------------------------------------------------------

    // small helper for approximate equality
    fn assert_almost_eq(a: f64, b: f64, tol: f64) {
        let diff = (a - b).abs();
        assert!(
            diff <= tol,
            "assert_almost_eq failed: left={a}, right={b}, diff={diff}, tol={tol}"
        );
    }

    // A direct KT zero-order coder for DNA (A,C,G,T) matching the same updates as CTW nodes.
    // Returns bits per base (bpb). Skips non-ACGT. No context.
    fn kt0_bits_per_base_dna(seq: &[u8]) -> f64 {
        #[inline]
        fn nuc_to_sym(b: u8) -> Option<usize> {
            match b {
                b'A' | b'a' => Some(0),
                b'C' | b'c' => Some(1),
                b'G' | b'g' => Some(2),
                b'T' | b't' => Some(3),
                _ => None,
            }
        }
        let m = 4usize;
        let mut counts = [0u32; 4];
        let mut n_eff = 0usize;
        let mut sum_log2 = 0.0;

        for &b in seq {
            if let Some(sym) = nuc_to_sym(b) {
                let c_s = counts[sym] as f64;
                let n = (counts[0] + counts[1] + counts[2] + counts[3]) as f64;
                // KT predictive factor: (c_s + 1/2)/(N + m/2), m=4 -> denom=N+2
                let num = c_s + 0.5;
                let den = n + (m as f64) / 2.0;
                sum_log2 += (num / den).ln() / std::f64::consts::LN_2;
                counts[sym] += 1;
                n_eff += 1;
            }
        }
        if n_eff == 0 {
            return 0.0;
        }
        let bits = -sum_log2;
        bits / (n_eff as f64)
    }

    #[test]
    fn test_ctw_depth0_equals_kt0() {
        // Mixed DNA including repeated patterns
        let s = b"ACGTACGTACGTGGGGCCCCAAAATTTTACGT";
        let kt0 = kt0_bits_per_base_dna(s);
        let ctw0 = super::ctw_bits_per_base_dna(s, 0);
        assert_almost_eq(ctw0, kt0, 1e-12);
    }

    #[test]
    fn test_ctw_context_improves_periodic() {
        // Highly periodic sequence → depth 1 should compress better than depth 0
        let s = b"ACACACACACACACACACACACACACACACAC";
        let d0 = super::ctw_bits_per_base_dna(s, 0);
        let d1 = super::ctw_bits_per_base_dna(s, 1);
        assert!(
            d1 <= d0,
            "depth 1 should not be worse than depth 0 (d1={d1}, d0={d0})"
        );

        // Depth 2 usually improves (or equals) further for strict period-2
        let d2 = super::ctw_bits_per_base_dna(s, 2);
        assert!(
            d2 <= d1 + 1e-12,
            "depth 2 should not be worse than depth 1 (d2={d2}, d1={d1})"
        );
    }

    #[test]
    fn test_ctw_skips_ns_depth0() {
        // In depth 0, skipping 'N' and flushing context is irrelevant;
        // result should match running on the sequence with N's removed.
        let with_ns = b"ACGTNNNNACGTNNAC";
        let no_ns: Vec<u8> = with_ns.iter().cloned().filter(|&b| b != b'N').collect();

        let ctw0_with = super::ctw_bits_per_base_dna(with_ns, 0);
        let ctw0_without = super::ctw_bits_per_base_dna(&no_ns, 0);
        assert_almost_eq(ctw0_with, ctw0_without, 1e-12);

        // Sanity: equals KT0 as well
        let kt0 = kt0_bits_per_base_dna(&no_ns);
        assert_almost_eq(ctw0_without, kt0, 1e-12);
    }

    #[test]
    fn test_ctw_empty_and_all_non_acgt() {
        let empty: &[u8] = b"";
        let only_ns: &[u8] = b"NNNNNNNN";
        let only_masked: &[u8] = b"nnnnxxxxNNNN";

        assert_eq!(super::ctw_bits_per_base_dna(empty, 6), 0.0);
        assert_eq!(super::ctw_bits_per_base_dna(only_ns, 6), 0.0);
        assert_eq!(super::ctw_bits_per_base_dna(only_masked, 6), 0.0);
    }

    #[test]
    fn test_ctw_reasonable_ranges() {
        // All the same base should be highly compressible (< 0.5 bpb)
        let same = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let bpb_same_d0 = super::ctw_bits_per_base_dna(same, 0);
        let bpb_same_d4 = super::ctw_bits_per_base_dna(same, 4);
        assert!(
            bpb_same_d0 < 0.5,
            "expected very low bits for mono-base, got {bpb_same_d0}"
        );
        assert!(
            bpb_same_d4 <= bpb_same_d0 + 1e-12,
            "context should not hurt for mono-base"
        );

        // Balanced i.i.d.-looking string should be around ~2 bits/base (slightly > because KT)
        let iid = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let bpb_iid_d0 = super::ctw_bits_per_base_dna(iid, 0);
        assert!(
            bpb_iid_d0 > 1.5 && bpb_iid_d0 < 2.2,
            "expected near-2 bits for balanced iid-like string, got {bpb_iid_d0}"
        );
    }
}
