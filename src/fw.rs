use itertools::Itertools;
use std::io::prelude::*;
use std::{fs::File, io::BufWriter};

use crate::kmer_maps::{self, KmerMap, WriteArray, WriteKmerValues};
use crate::kmeru8;
use crate::seq_statsu8;

use anyhow::Result;
use bio::io::fasta;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use std::sync::mpsc::channel;

pub fn fasta_windows(
    matches: &clap::ArgMatches,
    mut window_file_0: BufWriter<File>,
    mut window_file_1: BufWriter<File>,
    mut window_file_2: BufWriter<File>,
    mut window_file_3: BufWriter<File>,
    mut window_file_4: BufWriter<File>,
) -> Result<()> {
    // get matches
    let input_fasta = matches.value_of("fasta").unwrap();
    let output = matches.value_of("output").unwrap();
    let window_size: usize = matches.value_of_t("window_size").unwrap();
    let masked = matches.is_present("masked");
    let description = matches.is_present("description");

    // compute the 2-4mer kmer maps once only
    // hard code with false until I decide how to deal with
    // canonical kmers
    let kmer_maps = kmer_maps::generate_kmer_maps(false);

    // channel for collecting output
    let (sender, receiver) = channel();

    // iterate over fasta to get number of sequences
    // for the progress bar
    let mut nb_reads = 0;
    // read in the fasta from file
    let mut reader = fasta::Reader::from_file(input_fasta)
        .expect("[-]\tPath invalid.")
        .records();
    while let Some(Ok(_record)) = reader.next() {
        nb_reads += 1;
    }

    let progress_bar = ProgressBar::new(nb_reads);
    let pb_style = ProgressStyle::with_template(
        "[+]\tProcessing records: {bar:40.cyan/blue} {pos:>7}/{len:12}",
    )?
    .progress_chars(">>-");
    progress_bar.set_style(pb_style);

    // second reader for the computation
    eprintln!("[+]\tReading fasta from file");
    let reader = fasta::Reader::from_file(input_fasta).expect("[-]\tPath invalid.");
    reader
        .records()
        .par_bridge()
        .for_each_with(sender, |s, record| {
            let fasta_record = record.expect("[-]\tError during fasta record parsing.");

            // for the stats at the end.
            // initiate counters for the windows
            let mut start = 0;
            // https://github.com/tolkit/fasta_windows/issues/9
            // if the record length < window size, set end to the
            let mut end = match fasta_record.seq().len() < window_size {
                true => fasta_record.seq().len(),
                false => window_size,
            };

            // begin sliding windows
            // consider changing this to chunks_exact?
            let windows = fasta_record.seq().chunks(window_size);

            for win in windows {
                let seq_stats = seq_statsu8::seq_stats(win, masked);

                // unpack values
                let kmer_stats = kmeru8::kmer_diversity(win, kmer_maps.clone());

                // get description if present
                let desc = match fasta_record.desc() {
                    Some(d) => d.to_string(),
                    None => "No description.".to_string(),
                };

                s.send(Entry {
                    id: fasta_record.id().to_string(),
                    desc,
                    start,
                    end,
                    nuc_counts: seq_stats.nuc_counts,
                    gc_proportion: seq_stats.gc_proportion,
                    gc_skew: seq_stats.gc_skew,
                    at_skew: seq_stats.at_skew,
                    shannon_entropy: seq_stats.shannon_entropy,
                    g_s: seq_stats.g_s,
                    c_s: seq_stats.c_s,
                    a_s: seq_stats.a_s,
                    t_s: seq_stats.t_s,
                    n_s: seq_stats.n_s,
                    masked: seq_stats.masked,
                    cpg_s: ((*kmer_stats.di_freq.get(6).unwrap_or(&0) as f32) / seq_stats.len),
                    dinucleotides: kmer_stats.dinucleotides,
                    trinucleotides: kmer_stats.trinucleotides,
                    tetranucleotides: kmer_stats.tetranucleotides,
                    divalues: kmer_stats.di_freq,
                    trivalues: kmer_stats.tri_freq,
                    tetravalues: kmer_stats.tetra_freq,
                })
                .unwrap();

                // re-set the counter if counter > length of current sequence

                if end < fasta_record.seq().len() {
                    start += window_size;
                    end += window_size;

                    // now check if end overshoots the actual sequence length
                    // issue #8
                    if end > fasta_record.seq().len() {
                        end = fasta_record.seq().len()
                    }
                } else {
                    start = 0;
                    end = window_size;
                }
            }
            progress_bar.inc(1);
        });
    progress_bar.finish();
    let mut entries: Vec<Entry> = receiver.iter().collect();
    // parallel iteration messes up the order of scaffold ID's, so we fix that here.
    // however, within ID's windows should be ordered.
    entries.sort_by_key(|x| x.id.clone());

    let mut entry_writer = Output(entries);

    eprintln!("[+]\tWriting output to files");

    entry_writer.write_windows(&mut window_file_0, description)?;
    entry_writer.write_kmers(
        &mut window_file_1,
        &mut window_file_2,
        &mut window_file_3,
        &mut window_file_4,
        kmer_maps,
        description,
    )?;

    eprintln!("[+]\tOutput written to directory: ./fw_out/{}", output);

    Ok(())
}

// the output struct
pub struct Entry {
    // the id from the fasta file
    pub id: String,
    // the description from the fasta file
    pub desc: String,
    // the start of the window
    pub start: usize,
    // the end of the window
    pub end: usize,
    // the nucleotide counts
    pub nuc_counts: Vec<i32>,
    // the gc proportion
    pub gc_proportion: f32,
    // gc skew
    pub gc_skew: f32,
    // at skew
    pub at_skew: f32,
    // shannon entropy
    pub shannon_entropy: f64,
    // number of g's
    pub g_s: f32,
    // number of c's
    pub c_s: f32,
    // number of a's
    pub a_s: f32,
    // number of t's
    pub t_s: f32,
    // number of n's
    pub n_s: f32,
    // number of masked bases
    pub masked: f32,
    // number of cpg sites
    pub cpg_s: f32,
    // dinucleotide shannon entropy
    pub dinucleotides: f64,
    // trinucleotide shannon entropy
    pub trinucleotides: f64,
    // tetranucleotide shannon entropy
    pub tetranucleotides: f64,
    // the frequency distributions of 
    // each of the 3 kmer classes
    pub divalues: Vec<i32>,
    pub trivalues: Vec<i32>,
    pub tetravalues: Vec<i32>,
}

pub struct Output(Vec<Entry>);

impl Output {
    // write the windows file, optionally including a description
    pub fn write_windows(&mut self, file: &mut BufWriter<File>, description: bool) -> Result<()> {
        let header= match description {
            true => "ID\tdescription\tstart\tend\tGC_prop\tGC_skew\tAT_skew\tShannon_entropy\tProp_Gs\tProp_Cs\tProp_As\tProp_Ts\tProp_Ns\tProp_masked\tCpG_prop\tDinucleotide_Shannon\tTrinucleotide_Shannon\tTetranucleotide_Shannon".to_string(),
            false => "ID\tstart\tend\tGC_prop\tGC_skew\tAT_skew\tShannon_entropy\tProp_Gs\tProp_Cs\tProp_As\tProp_Ts\tProp_Ns\tProp_masked\tCpG_prop\tDinucleotide_Shannon\tTrinucleotide_Shannon\tTetranucleotide_Shannon".to_string()
        };

        writeln!(file, "{header}")?;

        for Entry {
            id,
            desc,
            start,
            end,
            nuc_counts: _,
            gc_proportion,
            gc_skew,
            at_skew,
            shannon_entropy,
            g_s,
            c_s,
            a_s,
            t_s,
            n_s,
            masked,
            cpg_s,
            dinucleotides,
            trinucleotides,
            tetranucleotides,
            divalues: _,
            trivalues: _,
            tetravalues: _,
        } in &self.0
        {
            let desc = match description {
                true => format!("{desc}\t"),
                false => String::new(),
            };
            writeln!(
                file,
                "{id}\t{desc}{start}\t{end}\t{gc_proportion:.3}\t{gc_skew:.3}\t{at_skew:.3}\t{shannon_entropy:.3}\t{g_s:.3}\t{c_s:.3}\t{a_s:.3}\t{t_s:.3}\t{n_s:.3}\t{masked:.3}\t{cpg_s:.3}\t{dinucleotides:.3}\t{trinucleotides:.3}\t{tetranucleotides:.3}",
            )?;
        }
        file.flush()?;

        Ok(())
    }

    pub fn write_kmers(
        &mut self,
        file1: &mut BufWriter<File>,
        file2: &mut BufWriter<File>,
        file3: &mut BufWriter<File>,
        file4: &mut BufWriter<File>,
        kmer_maps: Vec<KmerMap>,
        description: bool,
    ) -> Result<()> {
        // these are the arrays, tab separated bed-like format.
        // unpack kmer_maps
        match kmer_maps.as_slice() {
            [two, three, four] => {
                // headers for all
                let header = match description {
                    true => "ID\tdescription\tstart\tend\t".to_string(),
                    false => "ID\tstart\tend\t".to_string(),
                };

                // headers for mononucs
                writeln!(file1, "{header}A\tC\tG\tT\tN")?;
                // headers for dinucs
                let mut dinuc_headers = Vec::new();
                for key in two.map.keys().sorted() {
                    dinuc_headers.push(key)
                }
                let dinucs = WriteKmerValues(dinuc_headers);
                writeln!(file2, "{header}{dinucs}")?;
                // headers for trinucs
                let mut trinuc_headers = Vec::new();
                for key in three.map.keys().sorted() {
                    trinuc_headers.push(key)
                }
                let trinucs = WriteKmerValues(trinuc_headers);
                writeln!(file3, "{header}{trinucs}")?;
                // headers for tetranucs
                let mut tetranuc_headers = Vec::new();
                for key in four.map.keys().sorted() {
                    tetranuc_headers.push(key)
                }
                let tetranucs = WriteKmerValues(tetranuc_headers);
                writeln!(file4, "{header}{tetranucs}")?;

                for Entry {
                    id,
                    desc,
                    start,
                    end,
                    nuc_counts,
                    gc_proportion: _,
                    gc_skew: _,
                    at_skew: _,
                    shannon_entropy: _,
                    g_s: _,
                    c_s: _,
                    a_s: _,
                    t_s: _,
                    n_s: _,
                    masked: _,
                    cpg_s: _,
                    dinucleotides: _,
                    trinucleotides: _,
                    tetranucleotides: _,
                    divalues,
                    trivalues,
                    tetravalues,
                } in &self.0
                {
                    let desc = match description {
                        true => format!("{desc}\t"),
                        false => String::new(),
                    };

                    let nuc_counts = WriteArray(nuc_counts.clone());
                    writeln!(file1, "{id}\t{desc}{start}\t{end}\t{nuc_counts}")?;

                    let divalues_vec = WriteArray(divalues.clone());
                    writeln!(file2, "{id}\t{desc}{start}\t{end}\t{divalues_vec}",)?;

                    let trivalues_vec = WriteArray(trivalues.clone());
                    writeln!(file3, "{id}\t{desc}{start}\t{end}\t{trivalues_vec}",)?;

                    let tetravalues_vec = WriteArray(tetravalues.clone());
                    writeln!(file4, "{id}\t{desc}{start}\t{end}\t{tetravalues_vec}",)?;
                }

                file1.flush()?;
                file2.flush()?;
                file3.flush()?;
                file4.flush()?;
            }
            [..] => {} // Make the patterns exhaustive
        }
        Ok(())
    }
}
