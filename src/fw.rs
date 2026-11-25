use std::io::prelude::*;
use std::path::PathBuf;
use std::{fs::File, io::BufWriter};

use crate::kmeru8::{self, CPG_INDEX};
use crate::pipeline::{self, push_uint, Window};
use crate::seq_statsu8;

use anyhow::Result;

/// Context depth used for the CTW column.
pub const CTW_DEPTH: usize = 6;

pub fn fasta_windows(
    matches: &clap::ArgMatches,
    window_file_0: BufWriter<File>,
    window_file_1: BufWriter<File>,
    window_file_2: BufWriter<File>,
    window_file_3: BufWriter<File>,
    window_file_4: BufWriter<File>,
) -> Result<()> {
    // get matches
    let input_fasta = matches
        .get_one::<PathBuf>("fasta")
        .expect("handled by clap");
    let output = matches
        .get_one::<PathBuf>("output")
        .expect("handled by clap");
    let window_size = matches.get_one::<usize>("window_size").cloned().unwrap();
    let masked = matches.get_one::<bool>("masked").cloned().unwrap();
    let description = matches.get_one::<bool>("description").cloned().unwrap();
    let ctw = matches.get_one::<bool>("ctw").cloned().unwrap();

    let mut writers = vec![
        window_file_0,
        window_file_1,
        window_file_2,
        window_file_3,
        window_file_4,
    ];
    write_headers(&mut writers, description, ctw)?;

    eprintln!("[+]\tReading fasta from file");
    let mut writers = pipeline::run(input_fasta, window_size, writers, |window, bufs| {
        format_window(window, bufs, masked, description, ctw)
    })?;
    for w in &mut writers {
        w.flush()?;
    }

    eprintln!(
        "[+]\tOutput written to directory: ./fw_out/{}",
        output.display()
    );

    Ok(())
}

fn write_headers<W: Write>(writers: &mut [W], description: bool, ctw: bool) -> Result<()> {
    let header = match (description, ctw) {
        (true, true) => "ID\tdescription\tstart\tend\tGC_prop\tGC_skew\tAT_skew\tShannon_entropy\tctw\tProp_Gs\tProp_Cs\tProp_As\tProp_Ts\tProp_Ns\tProp_masked\tCpG_prop\tDinucleotide_Shannon\tTrinucleotide_Shannon\tTetranucleotide_Shannon",
        (false, true) => "ID\tstart\tend\tGC_prop\tGC_skew\tAT_skew\tShannon_entropy\tctw\tProp_Gs\tProp_Cs\tProp_As\tProp_Ts\tProp_Ns\tProp_masked\tCpG_prop\tDinucleotide_Shannon\tTrinucleotide_Shannon\tTetranucleotide_Shannon",
        (true, false) => "ID\tdescription\tstart\tend\tGC_prop\tGC_skew\tAT_skew\tShannon_entropy\tProp_Gs\tProp_Cs\tProp_As\tProp_Ts\tProp_Ns\tProp_masked\tCpG_prop\tDinucleotide_Shannon\tTrinucleotide_Shannon\tTetranucleotide_Shannon",
        (false, false) => "ID\tstart\tend\tGC_prop\tGC_skew\tAT_skew\tShannon_entropy\tProp_Gs\tProp_Cs\tProp_As\tProp_Ts\tProp_Ns\tProp_masked\tCpG_prop\tDinucleotide_Shannon\tTrinucleotide_Shannon\tTetranucleotide_Shannon",
    };
    writeln!(writers[0], "{header}")?;

    // the k-mer arrays: tab separated, bed-like format
    let prefix = match description {
        true => "ID\tdescription\tstart\tend\t",
        false => "ID\tstart\tend\t",
    };
    writeln!(writers[1], "{prefix}A\tC\tG\tT\tN")?;
    // writers[2..=4] hold the di-, tri- and tetranucleotide arrays
    for (w, k) in writers[2..].iter_mut().zip(2..=4) {
        writeln!(w, "{prefix}{}", kmeru8::kmer_names(k).join("\t"))?;
    }
    Ok(())
}

/// Append one window's rows to the five output buffers: windows, mono-, di-,
/// tri- and tetranucleotides.
fn format_window(w: &Window, bufs: &mut [Vec<u8>], masked: bool, description: bool, ctw: bool) {
    let stats = seq_statsu8::seq_stats(w.seq, masked);
    let kmers = kmeru8::count_kmers(w.seq);
    let cpg = kmers.di[CPG_INDEX] as f32 / stats.len;

    // "ID\t[description\t]start\tend" on every row
    for buf in bufs.iter_mut() {
        buf.extend_from_slice(w.record.id.as_bytes());
        buf.push(b'\t');
        if description {
            let desc = w.record.desc.as_deref().unwrap_or("No description.");
            buf.extend_from_slice(desc.as_bytes());
            buf.push(b'\t');
        }
        push_uint(buf, w.start as u64);
        buf.push(b'\t');
        push_uint(buf, w.end as u64);
    }

    // Vec<u8>'s io::Write never fails
    let b = &mut bufs[0];
    let _ = write!(
        b,
        "\t{:.3}\t{:.3}\t{:.3}\t{:.3}",
        stats.gc_proportion, stats.gc_skew, stats.at_skew, stats.shannon_entropy
    );
    if ctw {
        let _ = write!(
            b,
            "\t{:.3}",
            kmeru8::ctw_bits_per_base_dna(w.seq, CTW_DEPTH)
        );
    }
    let _ = writeln!(
        b,
        "\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}",
        stats.g_s,
        stats.c_s,
        stats.a_s,
        stats.t_s,
        stats.n_s,
        stats.masked,
        cpg,
        kmers.shannon_di(),
        kmers.shannon_tri(),
        kmers.shannon_tetra()
    );

    push_counts(&mut bufs[1], &stats.nuc_counts);
    push_counts(&mut bufs[2], &kmers.di);
    push_counts(&mut bufs[3], &kmers.tri);
    push_counts(&mut bufs[4], &kmers.tetra);
}

fn push_counts(buf: &mut Vec<u8>, counts: &[u64]) {
    for &c in counts {
        buf.push(b'\t');
        push_uint(buf, c);
    }
    buf.push(b'\n');
}
