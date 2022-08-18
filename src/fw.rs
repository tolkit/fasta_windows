use itertools::Itertools;
use std::io::prelude::*;
use std::{fs::File, io::BufWriter};

use crate::kmer_maps::{self, KmerMap, WriteArray, WriteKmerValues};
use crate::kmeru8;
use crate::seq_statsu8;

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
) -> std::io::Result<()> {
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
    let pb_style = ProgressStyle::default_bar()
        .template("[+]\tProcessing records: {bar:40.cyan/blue} {pos:>7}/{len:12}")
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
    pub id: String,
    pub desc: String,
    pub start: usize,
    pub end: usize,
    pub nuc_counts: Vec<i32>,
    pub gc_proportion: f32,
    pub gc_skew: f32,
    pub at_skew: f32,
    pub shannon_entropy: f64,
    pub g_s: f32,
    pub c_s: f32,
    pub a_s: f32,
    pub t_s: f32,
    pub n_s: f32,
    pub dinucleotides: f64,
    pub trinucleotides: f64,
    pub tetranucleotides: f64,
    pub divalues: Vec<i32>,
    pub trivalues: Vec<i32>,
    pub tetravalues: Vec<i32>,
}

pub struct Output(Vec<Entry>);

impl Output {
    // write the windows file, optionally including a description
    pub fn write_windows(
        &mut self,
        file: &mut BufWriter<File>,
        description: bool,
    ) -> std::io::Result<()> {
        let header;

        match description {
            true => header = "ID\tdescription\tstart\tend\tGC_prop\tGC_skew\tAT_skew\tShannon_entropy\tProp_Gs\tProp_Cs\tProp_As\tProp_Ts\tProp_Ns\tDinucleotide_Shannon\tTrinucleotide_Shannon\tTetranucleotide_Shannon".to_string(),
            false => header = "ID\tstart\tend\tGC_prop\tGC_skew\tAT_skew\tShannon_entropy\tProp_Gs\tProp_Cs\tProp_As\tProp_Ts\tProp_Ns\tDinucleotide_Shannon\tTrinucleotide_Shannon\tTetranucleotide_Shannon".to_string()
        }

        writeln!(file, "{}", header)
            .unwrap_or_else(|_| eprintln!("[-]\tError in writing to file."));

        for i in &self.0 {
            let desc = match description {
                true => format!("{}\t", i.desc),
                false => format!(""),
            };
            writeln!(
                file,
                "{}\t{}{}\t{}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}",
                i.id,
                desc,
                i.start,
                i.end,
                i.gc_proportion,
                i.gc_skew,
                i.at_skew,
                i.shannon_entropy,
                i.g_s,
                i.c_s,
                i.a_s,
                i.t_s,
                i.n_s,
                i.dinucleotides,
                i.trinucleotides,
                i.tetranucleotides,
            )
            .unwrap_or_else(|_| eprintln!("[-]\tError in writing to file."));
        }
        file.flush().unwrap();

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
    ) -> std::io::Result<()> {
        // these are the arrays, tab separated bed-like format.
        // unpack kmer_maps
        match kmer_maps.as_slice() {
            [two, three, four] => {
                // headers for all
                let header: String;
                match description {
                    true => header = "ID\tdescription\tstart\tend\t".to_string(),
                    false => header = "ID\tstart\tend\t".to_string(),
                }

                // headers for mononucs
                writeln!(file1, "{}A\tC\tG\tT", header)
                    .unwrap_or_else(|_| eprintln!("[-]\tError in writing to file."));
                // headers for dinucs
                let mut dinuc_headers = Vec::new();
                for key in two.map.keys().sorted() {
                    dinuc_headers.push(key)
                }
                writeln!(file2, "{}{}", header, WriteKmerValues(dinuc_headers))
                    .unwrap_or_else(|_| eprintln!("[-]\tError in writing to file."));
                // headers for trinucs
                let mut trinuc_headers = Vec::new();
                for key in three.map.keys().sorted() {
                    trinuc_headers.push(key)
                }
                writeln!(file3, "{}{}", header, WriteKmerValues(trinuc_headers))
                    .unwrap_or_else(|_| eprintln!("[-]\tError in writing to file."));
                // headers for tetranucs
                let mut tetranuc_headers = Vec::new();
                for key in four.map.keys().sorted() {
                    tetranuc_headers.push(key)
                }
                writeln!(file4, "{}{}", header, WriteKmerValues(tetranuc_headers))
                    .unwrap_or_else(|_| eprintln!("[-]\tError in writing to file."));

                for i in &self.0 {
                    let desc = match description {
                        true => format!("{}\t", i.desc),
                        false => format!(""),
                    };

                    writeln!(
                        file1,
                        "{}\t{}{}\t{}\t{}",
                        i.id,
                        desc,
                        i.start,
                        i.end,
                        WriteArray(i.nuc_counts.clone())
                    )
                    .unwrap_or_else(|_| eprintln!("[-]\tError in writing to file."));

                    writeln!(
                        file2,
                        "{}\t{}{}\t{}\t{}",
                        i.id,
                        desc,
                        i.start,
                        i.end,
                        WriteArray(i.divalues.clone())
                    )
                    .unwrap_or_else(|_| eprintln!("[-]\tError in writing to file."));

                    writeln!(
                        file3,
                        "{}\t{}{}\t{}\t{}",
                        i.id,
                        desc,
                        i.start,
                        i.end,
                        WriteArray(i.trivalues.clone())
                    )
                    .unwrap_or_else(|_| eprintln!("[-]\tError in writing to file."));

                    writeln!(
                        file4,
                        "{}\t{}{}\t{}\t{}",
                        i.id,
                        desc,
                        i.start,
                        i.end,
                        WriteArray(i.tetravalues.clone())
                    )
                    .unwrap_or_else(|_| eprintln!("[-]\tError in writing to file."));
                }

                file1.flush().unwrap();
                file2.flush().unwrap();
                file3.flush().unwrap();
                file4.flush().unwrap();
            }
            [..] => {} // Make the patterns exhaustive
        }
        Ok(())
    }
}
