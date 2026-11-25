use std::io::prelude::*;
use std::path::PathBuf;
use std::{fs::File, io::BufWriter};

use anyhow::Result;

use crate::fw::CTW_DEPTH;
use crate::pipeline::{self, push_uint};

// Bins each byte into one of 6 classes: A=0 C=1 G=2 T=3 N=4 other=5.
// Lowercase acgtn are folded to their uppercase equivalents.
const NUC_LUT: [u8; 256] = build_nuc_lut();
const fn build_nuc_lut() -> [u8; 256] {
    let mut lut = [5u8; 256];
    lut[b'A' as usize] = 0;
    lut[b'a' as usize] = 0;
    lut[b'C' as usize] = 1;
    lut[b'c' as usize] = 1;
    lut[b'G' as usize] = 2;
    lut[b'g' as usize] = 2;
    lut[b'T' as usize] = 3;
    lut[b't' as usize] = 3;
    lut[b'N' as usize] = 4;
    lut[b'n' as usize] = 4;
    lut
}

// Masked mode: only uppercase ACGTN are counted; everything else goes to bin 6 (skipped).
const MASKED_LUT: [u8; 256] = build_masked_lut();
const fn build_masked_lut() -> [u8; 256] {
    let mut lut = [6u8; 256];
    lut[b'A' as usize] = 0;
    lut[b'C' as usize] = 1;
    lut[b'G' as usize] = 2;
    lut[b'T' as usize] = 3;
    lut[b'N' as usize] = 4;
    lut
}

/// Shannon entropy of a DNA window using a compile-time LUT.
///
/// Counts 6 bins (ACGTN + other); ambiguous IUPAC codes (W, S, M, …) are
/// collapsed into "other". Four interleaved count tables avoid store-forwarding
/// stalls on runs of the same base.
///
/// masked=false: lowercase bases fold to uppercase and are counted normally.
/// masked=true:  only uppercase ACGTN contribute; denominator is their count.
#[inline(always)]
pub fn entropy_fast(dna: &[u8], masked: bool) -> f64 {
    let lut = if masked { &MASKED_LUT } else { &NUC_LUT };
    let mut tables = [[0u64; 8]; 4];
    let mut quads = dna.chunks_exact(4);
    for q in &mut quads {
        tables[0][lut[q[0] as usize] as usize] += 1;
        tables[1][lut[q[1] as usize] as usize] += 1;
        tables[2][lut[q[2] as usize] as usize] += 1;
        tables[3][lut[q[3] as usize] as usize] += 1;
    }
    for &b in quads.remainder() {
        tables[0][lut[b as usize] as usize] += 1;
    }
    let mut counts = [0u64; 6];
    for (i, c) in counts.iter_mut().enumerate() {
        *c = tables[0][i] + tables[1][i] + tables[2][i] + tables[3][i];
    }
    let total: u64 = counts.iter().sum();
    if total == 0 {
        return 0.0;
    }
    let n = total as f64;
    let mut entropy = 0f64;
    for &c in counts.iter().filter(|&&c| c > 0) {
        let p = c as f64 / n;
        entropy -= p * p.log2();
    }
    entropy
}

/// Fast path that only computes Shannon entropy and CTW and writes a BED file.
///
/// Uses the same streaming pipeline as the main mode: records are read in
/// batches, windows are processed in parallel across and within records, and
/// rows are written in FASTA order.
pub fn entropy_windows(matches: &clap::ArgMatches, bed_file: BufWriter<File>) -> Result<()> {
    let input_fasta = matches
        .get_one::<PathBuf>("fasta")
        .expect("handled by clap");
    let output = matches
        .get_one::<PathBuf>("output")
        .expect("handled by clap");
    let window_size = matches.get_one::<usize>("window_size").cloned().unwrap();
    let masked = matches.get_one::<bool>("masked").cloned().unwrap();

    eprintln!("[+]\tReading fasta (entropy mode)");
    let mut writers = pipeline::run(input_fasta, window_size, vec![bed_file], |w, bufs| {
        let entropy = entropy_fast(w.seq, masked);
        let ctw = crate::kmeru8::ctw_bits_per_base_dna(w.seq, CTW_DEPTH);
        let b = &mut bufs[0];
        b.extend_from_slice(w.record.id.as_bytes());
        b.push(b'\t');
        push_uint(b, w.start as u64);
        b.push(b'\t');
        push_uint(b, w.end as u64);
        // Vec<u8>'s io::Write never fails
        let _ = writeln!(b, "\t{entropy:.6}\t{ctw:.6}");
    })?;
    writers[0].flush()?;

    eprintln!(
        "[+]\tOutput written to: ./fw_out/{}_entropy.bed",
        output.display()
    );

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::entropy_fast;

    #[test]
    fn test_entropy_uniform() {
        // Equal ACGT → max entropy = 2.0 bits
        let dna = b"ACGTACGTACGTACGT";
        let h = entropy_fast(dna, false);
        assert!((h - 2.0).abs() < 1e-10, "expected 2.0, got {h}");
    }

    #[test]
    fn test_entropy_homopolymer() {
        let dna = b"AAAAAAAAAAAAAAAA";
        let h = entropy_fast(dna, false);
        assert!((h - 0.0).abs() < 1e-10, "expected 0.0, got {h}");
    }

    #[test]
    fn test_entropy_case_folding() {
        // Lowercase and uppercase should give the same result
        let upper = b"ACGTACGT";
        let lower = b"acgtacgt";
        let h_upper = entropy_fast(upper, false);
        let h_lower = entropy_fast(lower, false);
        assert!((h_upper - h_lower).abs() < 1e-10);
    }

    #[test]
    fn test_entropy_masked_excludes_lowercase() {
        // With masked=true, only uppercase ACGT count.
        // All-lowercase → total=0 → entropy=0.0
        let dna = b"acgtacgt";
        let h = entropy_fast(dna, true);
        assert_eq!(h, 0.0);
        // mixed: only the upper case half counts
        assert_eq!(entropy_fast(b"ACGTacgtRYK", true), 2.0);
    }

    #[test]
    fn test_entropy_odd_lengths() {
        // lengths that are not a multiple of the 4-way unroll
        assert_eq!(entropy_fast(b"AC", false), 1.0);
        assert_eq!(entropy_fast(b"ACGTA", false), entropy_fast(b"AACGT", false));
    }

    #[test]
    fn test_entropy_empty() {
        let h = entropy_fast(b"", false);
        assert_eq!(h, 0.0);
    }

    #[test]
    fn test_entropy_all_n() {
        // All-N → single class → entropy = 0
        let dna = b"NNNNNNNN";
        let h = entropy_fast(dna, false);
        assert_eq!(h, 0.0);
    }
}
