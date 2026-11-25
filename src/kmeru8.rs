use std::cell::RefCell;
use std::sync::OnceLock;

/// 2-bit code for each byte: A=0 C=1 G=2 T=3 (either case); 4 for anything else.
pub const NUC_2BIT: [u8; 256] = {
    let mut t = [4u8; 256];
    t[b'A' as usize] = 0;
    t[b'a' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'c' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'g' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b't' as usize] = 3;
    t
};

/// Index of the CG dinucleotide in [`KmerCounts::di`] (C=1, G=2 → 1*4 + 2).
pub const CPG_INDEX: usize = 6;

/// Di-, tri- and tetranucleotide counts of a window.
///
/// Index `i` of each array is the k-mer whose bases, read as base-4 digits
/// (A=0 C=1 G=2 T=3, first base most significant), spell `i`. That is
/// lexicographic order, AA..TT, matching the output column order.
pub struct KmerCounts {
    pub di: [u64; 16],
    pub tri: [u64; 64],
    pub tetra: [u64; 256],
}

impl KmerCounts {
    pub fn shannon_di(&self) -> f64 {
        shannon(&self.di)
    }
    pub fn shannon_tri(&self) -> f64 {
        shannon(&self.tri)
    }
    pub fn shannon_tetra(&self) -> f64 {
        shannon(&self.tetra)
    }
}

/// Count all 2-, 3- and 4-mers in `dna` in one pass.
///
/// Case is ignored. A k-mer is only counted if all k bases are A, C, G or T;
/// any other byte (N, IUPAC codes, ...) breaks the run.
pub fn count_kmers(dna: &[u8]) -> KmerCounts {
    let mut k = KmerCounts {
        di: [0; 16],
        tri: [0; 64],
        tetra: [0; 256],
    };
    let mut code = 0usize; // last four bases, most recent in the low bits
    let mut run = 0u32; // consecutive ACGT bases ending here
    for &b in dna {
        let sym = NUC_2BIT[b as usize];
        if sym > 3 {
            run = 0;
            continue;
        }
        code = ((code << 2) | sym as usize) & 0xFF;
        run += 1;
        if run >= 4 {
            k.di[code & 0xF] += 1;
            k.tri[code & 0x3F] += 1;
            k.tetra[code] += 1;
        } else if run == 3 {
            k.di[code & 0xF] += 1;
            k.tri[code & 0x3F] += 1;
        } else if run == 2 {
            k.di[code & 0xF] += 1;
        }
    }
    k
}

/// Header names for the k-mers in [`KmerCounts`] order.
pub fn kmer_names(k: usize) -> Vec<String> {
    (0..1usize << (2 * k))
        .map(|i| {
            (0..k)
                .rev()
                .map(|pos| ['A', 'C', 'G', 'T'][(i >> (2 * pos)) & 3])
                .collect()
        })
        .collect()
}

/// Shannon diversity (bits) of a frequency distribution.
pub fn shannon(counts: &[u64]) -> f64 {
    let total: u64 = counts.iter().sum();
    let mut diversity = 0f64;
    for &count in counts.iter().filter(|&&c| c > 0) {
        let p = count as f64 / total as f64;
        diversity -= p * p.log2();
    }
    diversity
}

// ---------------------------------------------------------------------------
// Context Tree Weighting
// ---------------------------------------------------------------------------

/// Largest supported CTW context depth (the tree is a flat array of
/// (4^(depth+1) - 1) / 3 nodes per thread: 5,461 at depth 6, 1.4M at depth 10).
pub const MAX_CTW_DEPTH: usize = 10;

/// Entries in the tabulated KT log-probabilities; larger counts use Stirling's series.
const KT_TABLE_LEN: usize = 1 << 16;

/// f[c] = log2 prod_{j<c} (j + 1/2) and g[n] = log2 prod_{i<n} (i + 2): the KT
/// code length of a node with counts c_a and total n is sum_a f[c_a] - g[n].
struct KtTables {
    f: Vec<f64>,
    g: Vec<f64>,
}

fn kt_tables() -> &'static KtTables {
    static TABLES: OnceLock<KtTables> = OnceLock::new();
    TABLES.get_or_init(|| {
        let mut f = vec![0.0; KT_TABLE_LEN];
        let mut g = vec![0.0; KT_TABLE_LEN];
        for i in 1..KT_TABLE_LEN {
            f[i] = f[i - 1] + ((i - 1) as f64 + 0.5).log2();
            g[i] = g[i - 1] + ((i - 1) as f64 + 2.0).log2();
        }
        KtTables { f, g }
    })
}

/// ln Γ(x) by Stirling's series; accurate to f64 precision for x >= KT_TABLE_LEN.
fn ln_gamma_large(x: f64) -> f64 {
    let inv = 1.0 / x;
    let inv2 = inv * inv;
    (x - 0.5) * x.ln() - x
        + 0.5 * (2.0 * std::f64::consts::PI).ln()
        + inv * (1.0 / 12.0 - inv2 * (1.0 / 360.0 - inv2 * (1.0 / 1260.0 - inv2 / 1680.0)))
}

impl KtTables {
    #[inline(always)]
    fn f(&self, c: u32) -> f64 {
        match self.f.get(c as usize) {
            Some(&v) => v,
            // log2(Γ(c + 1/2) / Γ(1/2)), with ln Γ(1/2) = ln(π)/2
            None => {
                (ln_gamma_large(c as f64 + 0.5) - 0.5 * std::f64::consts::PI.ln())
                    * std::f64::consts::LOG2_E
            }
        }
    }

    #[inline(always)]
    fn g(&self, n: u32) -> f64 {
        match self.g.get(n as usize) {
            Some(&v) => v,
            // log2(Γ(n + 2) / Γ(2))
            None => ln_gamma_large(n as f64 + 2.0) * std::f64::consts::LOG2_E,
        }
    }

    #[inline(always)]
    fn kt(&self, c: [u32; 4]) -> f64 {
        self.f(c[0]) + self.f(c[1]) + self.f(c[2]) + self.f(c[3])
            - self.g(c[0] + c[1] + c[2] + c[3])
    }
}

/// log2(2^a / 2 + 2^b / 2), the CTW mixture with beta = 1/2.
#[inline(always)]
fn log2_mix_half(a: f64, b: f64) -> f64 {
    let m = a.max(b);
    let d = (a - b).abs();
    if d > 50.0 {
        m - 1.0
    } else {
        m - 1.0 + (-d).exp2().ln_1p() * std::f64::consts::LOG2_E
    }
}

/// Per-thread CTW working memory, reused across windows.
struct CtwScratch {
    depth: usize,
    /// KT counts per node. Node 0 is the root; the children of node i are
    /// 4i+1 ..= 4i+4, so every child has a larger index than its parent.
    counts: Vec<[u32; 4]>,
    /// whether the node's most recent update was as the deepest node on its path
    leaf: Vec<bool>,
    /// log2 weighted probability
    log_w: Vec<f64>,
    /// nodes touched in the current window, by depth
    touched: Vec<Vec<u32>>,
}

impl CtwScratch {
    fn with_depth(depth: usize) -> Self {
        let nodes = ((1usize << (2 * (depth + 1))) - 1) / 3;
        CtwScratch {
            depth,
            counts: vec![[0; 4]; nodes],
            leaf: vec![false; nodes],
            log_w: vec![0.0; nodes],
            touched: vec![Vec::new(); depth + 1],
        }
    }

    fn bits_per_base(&mut self, dna: &[u8]) -> f64 {
        let depth = self.depth;
        let mask = (1usize << (2 * depth)) - 1;
        let mut context = 0usize; // previous `depth` symbols, most recent in the low bits
        let mut context_len = 0usize;
        let mut n_eff = 0usize;

        // Pass 1: per base, only add the symbol to the KT counts of each node
        // on its context path.
        for &b in dna {
            let sym = NUC_2BIT[b as usize];
            if sym > 3 {
                context_len = 0; // flush the context on N, IUPAC codes, etc.
                continue;
            }
            let sym = sym as usize;
            let mut node = 0usize;
            for d in 0..=context_len {
                if d > 0 {
                    node = 4 * node + 1 + ((context >> (2 * (d - 1))) & 3);
                }
                let c = &mut self.counts[node];
                if *c == [0; 4] {
                    self.touched[d].push(node as u32);
                }
                c[sym] += 1;
                self.leaf[node] = d == context_len;
            }
            n_eff += 1;
            context = ((context << 2) | sym) & mask;
            if context_len < depth {
                context_len += 1;
            }
        }
        if n_eff == 0 {
            return 0.0;
        }

        // Pass 2: weighted probabilities, bottom-up, once per touched node.
        let kt = kt_tables();
        for d in (0..=depth).rev() {
            for &node in &self.touched[d] {
                let node = node as usize;
                let log_kt = kt.kt(self.counts[node]);
                self.log_w[node] = if d == depth || self.leaf[node] {
                    log_kt
                } else {
                    let c = 4 * node + 1;
                    // untouched children have log_w == 0 (probability 1)
                    let children = 0.0
                        + self.log_w[c]
                        + self.log_w[c + 1]
                        + self.log_w[c + 2]
                        + self.log_w[c + 3];
                    log2_mix_half(log_kt, children)
                };
            }
        }
        let bits = -self.log_w[0];

        for d in 0..=depth {
            for &node in &self.touched[d] {
                self.counts[node as usize] = [0; 4];
                self.log_w[node as usize] = 0.0;
            }
            self.touched[d].clear();
        }
        bits / n_eff as f64
    }
}

thread_local! {
    static CTW_SCRATCH: RefCell<Option<CtwScratch>> = const { RefCell::new(None) };
}

/// Context Tree Weighting (CTW) code length for DNA (A,C,G,T), in bits per
/// effective base, using context depth `max_depth` (at most [`MAX_CTW_DEPTH`]).
/// Non-ACGT symbols are skipped and flush the context. With `max_depth == 0`
/// this is the zero-order KT code length.
///
/// - Willems, Shtarkov & Tjalkens (1995): CTW basics.
/// - Multinomial KT estimator with 1/2 pseudo-counts; mixture weight beta = 1/2.
/// - A node reached with an exhausted (flushed) context is updated as a leaf:
///   its weighted probability becomes its KT probability.
///
/// Rather than updating every node's weighted probability after every base,
/// the window is processed in two passes. A node's weighted probability only
/// changes when the node itself is updated, every update of a node also
/// updates all its ancestors (after it), and the root's per-base changes sum
/// to its final value. So the window's code length is the root's final
/// weighted probability, which depends only on each node's final KT counts and
/// on whether its last update was as a leaf. Pass 1 collects those; pass 2
/// evaluates the tree once, bottom-up. This needs one mixture per touched
/// node per window instead of one per node per base.
pub fn ctw_bits_per_base_dna(dna: &[u8], max_depth: usize) -> f64 {
    assert!(
        max_depth <= MAX_CTW_DEPTH,
        "CTW depth {max_depth} exceeds the maximum of {MAX_CTW_DEPTH}"
    );
    CTW_SCRATCH.with(|cell| {
        let mut scratch = cell.borrow_mut();
        if scratch.as_ref().is_none_or(|s| s.depth != max_depth) {
            *scratch = Some(CtwScratch::with_depth(max_depth));
        }
        scratch.as_mut().unwrap().bits_per_base(dna)
    })
}

pub fn reverse_complement(dna: &[u8]) -> Vec<u8> {
    dna.iter().rev().map(|&b| switch_base(b)).collect()
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
    use std::collections::HashMap;

    #[test]
    fn test_reverse_complement() {
        let short_dna_string = "AACCTTGG".as_bytes();

        let revcomp_dna = reverse_complement(short_dna_string);

        let should_be = "CCAAGGTT".as_bytes().to_vec();

        assert_eq!(should_be, revcomp_dna);
    }

    #[test]
    fn test_shannon_entropy() {
        // uniform over 16 dinucleotides: 4 bits
        assert_eq!(shannon(&[3; 16]), 4.0);
        // a single class: 0 bits
        assert_eq!(shannon(&[0, 7, 0, 0]), 0.0);
        // two equal classes: 1 bit
        assert_eq!(shannon(&[5, 0, 5, 0]), 1.0);
        // nothing counted
        assert_eq!(shannon(&[0; 16]), 0.0);
        let k = count_kmers(b"ACGTACGTACGTACGTACGT");
        // 4 distinct dinucleotides (AC, CG, GT, TA) in near-equal numbers
        assert!((k.shannon_di() - 2.0).abs() < 0.01);
    }

    #[test]
    fn test_kmer_names() {
        assert_eq!(kmer_names(2)[..5], ["AA", "AC", "AG", "AT", "CA"]);
        assert_eq!(kmer_names(2)[CPG_INDEX], "CG");
        assert_eq!(kmer_names(3).len(), 64);
        assert_eq!(kmer_names(4)[255], "TTTT");
    }

    // Deterministic pseudo-random sequences for the oracle tests.
    fn random_seq(seed: u64, len: usize, alphabet: &[u8]) -> Vec<u8> {
        let mut state = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (0..len)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                alphabet[(state % alphabet.len() as u64) as usize]
            })
            .collect()
    }

    /// Naive k-mer counting: upper-case each k-mer, skip any containing a
    /// non-ACGT byte.
    fn reference_kmers(dna: &[u8], k: usize) -> Vec<u64> {
        let names = kmer_names(k);
        let mut map: HashMap<Vec<u8>, u64> =
            names.iter().map(|n| (n.as_bytes().to_vec(), 0)).collect();
        for kmer in dna.windows(k) {
            let upper = kmer.to_ascii_uppercase();
            if let Some(c) = map.get_mut(&upper) {
                *c += 1;
            }
        }
        names.iter().map(|n| map[n.as_bytes()]).collect()
    }

    #[test]
    fn test_count_kmers_matches_reference() {
        for (seed, alphabet) in [
            (1, &b"ACGT"[..]),
            (2, b"ACGTacgt"),
            (3, b"ACGTNacgtn"),
            (4, b"AACCGGTTNRYKMSWacgt"),
            (5, b"AAAAAAAAAC"),
        ] {
            for len in [0, 1, 2, 3, 4, 5, 17, 1000] {
                let dna = random_seq(seed, len, alphabet);
                let k = count_kmers(&dna);
                assert_eq!(k.di.to_vec(), reference_kmers(&dna, 2));
                assert_eq!(k.tri.to_vec(), reference_kmers(&dna, 3));
                assert_eq!(k.tetra.to_vec(), reference_kmers(&dna, 4));
            }
        }
    }

    #[test]
    fn test_iupac_breaks_kmers() {
        let k = count_kmers(b"ACRGTyAC");
        let names = kmer_names(2);
        let seen: Vec<_> = (0..16)
            .filter(|&i| k.di[i] > 0)
            .map(|i| names[i].as_str())
            .collect();
        assert_eq!(seen, ["AC", "GT"]);
        assert_eq!(k.di[0b0001], 2); // AC twice
        assert_eq!(k.tri.iter().sum::<u64>(), 0);
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

    /// The original per-base, pointer-based CTW (fasta_windows 0.2.6), kept as
    /// the oracle for the two-pass implementation.
    fn reference_ctw(dna: &[u8], max_depth: usize) -> f64 {
        if max_depth == 0 {
            return kt0_bits_per_base_dna(dna);
        }
        fn nuc_to_sym(b: u8) -> Option<usize> {
            match b {
                b'A' | b'a' => Some(0),
                b'C' | b'c' => Some(1),
                b'G' | b'g' => Some(2),
                b'T' | b't' => Some(3),
                _ => None,
            }
        }
        fn log2(x: f64) -> f64 {
            x.ln() / std::f64::consts::LN_2
        }
        fn log2_sum_weighted(a_log2: f64, b_log2: f64, beta: f64) -> f64 {
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
            counts: [u32; 4],
            total: u32,
            log_p_kt: f64,
            log_w: f64,
            children: [Option<Box<Node>>; 4],
        }
        impl Node {
            fn new() -> Self {
                Node {
                    counts: [0; 4],
                    total: 0,
                    log_p_kt: 0.0,
                    log_w: 0.0,
                    children: [None, None, None, None],
                }
            }
            fn kt_update(&mut self, sym: usize) {
                let num = self.counts[sym] as f64 + 0.5;
                let den = self.total as f64 + 2.0;
                self.log_p_kt += (num / den).ln() / std::f64::consts::LN_2;
                self.counts[sym] += 1;
                self.total += 1;
            }
            fn update(&mut self, ctx: &[usize], sym: usize) {
                if let Some((&a, rest)) = ctx.split_first() {
                    self.children[a]
                        .get_or_insert_with(|| Box::new(Node::new()))
                        .update(rest, sym);
                    self.kt_update(sym);
                    let mut sum = 0.0;
                    for child in self.children.iter().flatten() {
                        sum += child.log_w;
                    }
                    self.log_w = log2_sum_weighted(self.log_p_kt, sum, 0.5);
                } else {
                    self.kt_update(sym);
                    self.log_w = self.log_p_kt;
                }
            }
        }
        let mut root = Node::new();
        let mut ctx: Vec<usize> = Vec::with_capacity(max_depth);
        let mut total = 0.0;
        let mut n_eff = 0usize;
        for &b in dna {
            let Some(sym) = nuc_to_sym(b) else {
                ctx.clear();
                continue;
            };
            let before = root.log_w;
            root.update(&ctx, sym);
            total += root.log_w - before;
            n_eff += 1;
            if ctx.len() == max_depth {
                ctx.pop();
            }
            ctx.insert(0, sym);
        }
        if n_eff == 0 {
            return 0.0;
        }
        -total / n_eff as f64
    }

    fn ctw_cases() -> Vec<Vec<u8>> {
        let mut cases = vec![
            vec![b'A'; 1000],
            b"ACGTTGA".repeat(140),
            b"ACACACAC".repeat(125),
            b"ACGTACGTACGTGGGGCCCCAAAATTTTACGT".to_vec(),
            b"ACGTNNNNACGTNNAC".to_vec(),
            b"nnnnxxxxNNNN".to_vec(),
            b"A".to_vec(),
            Vec::new(),
            random_seq(10, 1000, b"ACGT"),
            random_seq(11, 1000, b"ACGTacgtN"),
            random_seq(12, 300, b"AAAAAAAAACGTNR"),
        ];
        let mut v = b"AGCTTAGC".repeat(60);
        v.extend(random_seq(13, 400, b"ACGT"));
        cases.push(v);
        let mut v = random_seq(14, 500, b"ACGT");
        v.extend(b"TTAGGG".repeat(80));
        v.extend(random_seq(15, 20, b"ACGT"));
        cases.push(v);
        // frequent context flushes, so many inner nodes are updated as leaves
        let mut v = Vec::new();
        for i in 0..40 {
            v.extend(random_seq(16 + i, 17 + i as usize, b"ACGT"));
            v.push(if i % 3 == 0 { b'N' } else { b'n' });
        }
        cases.push(v);
        let mut v = b"AT".repeat(2000);
        v.extend(random_seq(60, 1000, b"ACGT"));
        cases.push(v);
        cases
    }

    #[test]
    fn test_ctw_matches_reference() {
        for depth in 0..=6 {
            for case in ctw_cases() {
                let expected = reference_ctw(&case, depth);
                let got = ctw_bits_per_base_dna(&case, depth);
                assert_almost_eq(got, expected, 1e-12);
            }
        }
    }

    #[test]
    fn test_ctw_scratch_reuse() {
        // the per-thread scratch must be fully reset between windows and depths
        let a = random_seq(100, 2000, b"ACGTN");
        let b = b"ACGTTGA".repeat(50);
        let first = ctw_bits_per_base_dna(&a, 6);
        let _ = ctw_bits_per_base_dna(&b, 6);
        let _ = ctw_bits_per_base_dna(&b, 3);
        assert_eq!(ctw_bits_per_base_dna(&a, 6).to_bits(), first.to_bits());
    }

    #[test]
    fn test_ctw_counts_beyond_kt_table() {
        // more symbols in a node than KT_TABLE_LEN, so Stirling's series is used
        let mut dna = vec![b'A'; KT_TABLE_LEN + 5000];
        dna.extend(random_seq(200, 20000, b"ACGT"));
        for depth in [0, 2] {
            assert_almost_eq(
                ctw_bits_per_base_dna(&dna, depth),
                reference_ctw(&dna, depth),
                1e-9,
            );
        }
        // the tabulated sums and Stirling's series agree where they meet
        let tables = kt_tables();
        let c = KT_TABLE_LEN as u32 - 1;
        let stirling_f = (ln_gamma_large(c as f64 + 0.5) - 0.5 * std::f64::consts::PI.ln())
            * std::f64::consts::LOG2_E;
        let stirling_g = ln_gamma_large(c as f64 + 2.0) * std::f64::consts::LOG2_E;
        assert_almost_eq(tables.f(c), stirling_f, 1e-6);
        assert_almost_eq(tables.g(c), stirling_g, 1e-6);
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

        // A balanced string should cost about 2 bits/base, plus the KT learning
        // cost of roughly (3/2)·log2(n)/n ≈ 0.23 bits/base at n = 32.
        let iid = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let bpb_iid_d0 = super::ctw_bits_per_base_dna(iid, 0);
        assert!(
            bpb_iid_d0 > 1.5 && bpb_iid_d0 < 2.3,
            "expected near-2 bits for balanced iid-like string, got {bpb_iid_d0}"
        );
    }
}
