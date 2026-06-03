use std::io::prelude::*;
use std::path::PathBuf;
use std::{fs::File, io::BufWriter};

use anyhow::Result;
use indicatif::{ProgressBar, ProgressStyle};
use needletail::parse_fastx_file;
use rayon::prelude::*;

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

// Masked mode: only uppercase ACGTN are counted; everything else maps to 255 (skip).
const MASKED_LUT: [u8; 256] = build_masked_lut();
const fn build_masked_lut() -> [u8; 256] {
    let mut lut = [255u8; 256];
    lut[b'A' as usize] = 0;
    lut[b'C' as usize] = 1;
    lut[b'G' as usize] = 2;
    lut[b'T' as usize] = 3;
    lut[b'N' as usize] = 4;
    lut
}

/// Shannon entropy of a DNA window using a compile-time LUT.
///
/// Uses 6-bin counting (ACGTN + other) over a single array that fits in one
/// cache line, enabling better auto-vectorisation than the 256-entry fallback.
/// Ambiguous IUPAC codes (W, S, M, …) are collapsed into the "other" bin.
///
/// masked=false: lowercase bases fold to uppercase and are counted normally.
/// masked=true:  only uppercase ACGTN contribute; denominator is their count.
#[inline(always)]
pub fn entropy_fast(dna: &[u8], masked: bool) -> f64 {
    let mut counts = [0u32; 6];
    if masked {
        for &byte in dna {
            let idx = MASKED_LUT[byte as usize];
            if idx < 255 {
                counts[idx as usize] += 1;
            }
        }
    } else {
        for &byte in dna {
            counts[NUC_LUT[byte as usize] as usize] += 1;
        }
    }
    let total: u32 = counts.iter().sum();
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

/// Fast path that only computes Shannon entropy and writes a BED file.
///
/// Architecture:
///   1. Read all records sequentially into memory (bio strips newlines for free).
///   2. Compute entropy in parallel with rayon's indexed par_iter + par_chunks:
///      - par_iter across chromosomes keeps all cores busy for many-scaffold genomes.
///      - par_chunks within each chromosome fully utilises cores for few-scaffold genomes.
///      - rayon's IndexedParallelIterator collect() preserves FASTA file order,
///        so no sort or MPSC channel is needed.
///   3. Write in order directly.
pub fn entropy_windows(matches: &clap::ArgMatches, mut bed_file: BufWriter<File>) -> Result<()> {
    let input_fasta = matches
        .get_one::<PathBuf>("fasta")
        .expect("handled by clap");
    let output = matches
        .get_one::<PathBuf>("output")
        .expect("handled by clap");
    let window_size = matches.get_one::<usize>("window_size").cloned().unwrap();
    let masked = matches.get_one::<bool>("masked").cloned().unwrap();

    let spinner = ProgressBar::new_spinner();
    spinner.set_style(ProgressStyle::with_template(
        "[+]\tReading: {spinner:.cyan} {pos} sequences",
    )?);

    eprintln!("[+]\tReading fasta (entropy mode)");
    // needletail's seq() strips newlines with SIMD memchr2 and returns Cow::Owned
    // for multi-line sequences (the common case), so into_owned() is a move not a copy.
    let mut records: Vec<(String, Vec<u8>)> = Vec::new();
    let mut reader = parse_fastx_file(input_fasta)?;
    while let Some(result) = reader.next() {
        let rec = result?;
        let raw_id = rec.id();
        let id_end = raw_id
            .iter()
            .position(|&b| b == b' ' || b == b'\t')
            .unwrap_or(raw_id.len());
        let id = String::from_utf8_lossy(&raw_id[..id_end]).into_owned();
        let seq = rec.seq().into_owned();
        records.push((id, seq));
        spinner.inc(1);
    }
    spinner.finish();

    // par_iter on a Vec is an IndexedParallelIterator: collect() preserves order.
    // par_chunks within each record gives intra-chromosome parallelism so a
    // genome dominated by a handful of large scaffolds still uses all cores.
    let results: Vec<Vec<(usize, usize, f64, f64)>> = records
        .par_iter()
        .map(|(_, seq)| {
            seq.par_chunks(window_size)
                .enumerate()
                .map(|(i, win)| {
                    let start = i * window_size;
                    let entropy = entropy_fast(win, masked);
                    let ctw = crate::kmeru8::ctw_bits_per_base_dna(win, 6);
                    (start, start + win.len(), entropy, ctw)
                })
                .collect()
        })
        .collect();

    eprintln!("[+]\tWriting BED output");
    for ((id, _), windows) in records.iter().zip(results.iter()) {
        for &(start, end, entropy, ctw) in windows {
            writeln!(
                bed_file,
                "{}\t{}\t{}\t{:.6}\t{:.6}",
                id, start, end, entropy, ctw
            )?;
        }
    }
    bed_file.flush()?;

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
