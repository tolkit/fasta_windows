// external imports
// std
use std::fs::{create_dir_all, File};
use std::io::prelude::*;
use std::io::BufWriter;
use std::sync::mpsc::channel;

// non-std
use bio::io::fasta;
use clap::{value_t, App, Arg};
use indicatif::{ProgressBar, ProgressStyle};
use itertools::Itertools;
use rayon::prelude::*;

// internal imports
use fasta_windows::kmer_maps::kmer_maps;
use fasta_windows::kmeru8::kmeru8;
use fasta_windows::seq_statsu8::seq_statsu8;

fn main() {
    // command line options
    let matches = App::new("Fasta windows")
        .version(clap::crate_version!())
        .author("Max Brown <mb39@sanger.ac.uk>")
        .about("Quickly compute statistics over a fasta file in windows.")
        .arg(
            Arg::with_name("fasta")
                .short("f")
                .long("fasta")
                .takes_value(true)
                .required(true)
                .help("The input fasta file."),
        )
        .arg(
            Arg::with_name("window_size")
                .short("w")
                .long("window_size")
                .help("Integer size of window for statistics to be computed over.")
                .takes_value(true)
                .default_value("1000"),
        )
        .arg(
            Arg::with_name("canonical_kmers")
                .short("c")
                .long("canonical_kmers")
                .help("Should the canonical kmers be calculated?"),
        )
        .arg(
            Arg::with_name("output")
                .short("o")
                .long("output")
                .help("Output filename for the CSV (without extension).")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("masked")
                .short("m")
                .long("masked")
                .help("Consider only uppercase nucleotides in the calculations."),
        )
        .get_matches();
    // parse command line options
    let input_fasta = matches.value_of("fasta").unwrap();
    let output = matches.value_of("output").unwrap();
    let window_size = value_t!(matches.value_of("window_size"), usize).unwrap_or_else(|e| e.exit());
    let canonical_kmers = matches.is_present("canonical_kmers");
    let masked = matches.is_present("masked");

    // create directory for output
    if let Err(e) = create_dir_all("./fw_out/") {
        eprintln!("[-]\tCreate directory error: {}", e.to_string());
    }

    // initiate the output TSV for windows over genome
    let output_file_1 = format!("./fw_out/{}{}", output, "_windows.tsv");
    let window_file = File::create(&output_file_1).unwrap();
    let mut window_file = BufWriter::new(window_file);

    writeln!(
        window_file,
        "ID\tstart\tend\tGC_prop\tGC_skew\tShannon_entropy\tProp_Gs\tProp_Cs\tProp_As\tProp_Ts\tProp_Ns\tDinucleotide_Shannon{arg}\tTrinucleotide_Shannon{arg}\tTetranucleotide_Shannon{arg}",
        arg = format!("_{}", canonical_kmers)
    )
    .unwrap();

    let output_file_2 = format!("./fw_out/{}{}", output, "_dinuc_windows.tsv");
    let window_file_2 = File::create(&output_file_2).unwrap();
    let mut window_file_2 = BufWriter::new(window_file_2);

    let output_file_3 = format!("./fw_out/{}{}", output, "_trinuc_windows.tsv");
    let window_file_3 = File::create(&output_file_3).unwrap();
    let mut window_file_3 = BufWriter::new(window_file_3);

    let output_file_4 = format!("./fw_out/{}{}", output, "_tetranuc_windows.tsv");
    let window_file_4 = File::create(&output_file_4).unwrap();
    let mut window_file_4 = BufWriter::new(window_file_4);

    // the output struct
    struct Output {
        id: String,
        start: usize,
        end: usize,
        gc_proportion: f32,
        gc_skew: f32,
        shannon_entropy: f64,
        g_s: f32,
        c_s: f32,
        a_s: f32,
        t_s: f32,
        n_s: f32,
        dinucleotides: f64,
        trinucleotides: f64,
        tetranucleotides: f64,
        divalues: Vec<i32>,
        trivalues: Vec<i32>,
        tetravalues: Vec<i32>,
    }

    // compute the 2-4mer kmer maps once only
    let kmer_maps = kmer_maps::generate_kmer_maps(canonical_kmers);

    // channel for collecting output
    let (sender, receiver) = channel();

    // iterate over fasta to get number of sequences
    // for the progress bar
    // this may waste a little time, but is there another option?
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
            // read_lengths.push(record.seq().len());
            // initiate a counter for the windows
            let mut counter = window_size;
            // begin sliding windows
            let windows = fasta_record.seq().chunks(window_size);

            for win in windows {
                let seq_stats = seq_statsu8::seq_stats(win, masked);

                // unpack values
                let kmer_stats = kmeru8::kmer_diversity(win, kmer_maps.clone(), canonical_kmers);

                s.send(Output {
                    id: fasta_record.id().to_string(),
                    start: counter - window_size,
                    end: counter,
                    gc_proportion: seq_stats.gc_proportion,
                    gc_skew: seq_stats.gc_skew,
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
                if counter < fasta_record.seq().len() {
                    counter += window_size
                } else {
                    counter = 0
                }
            }
            progress_bar.inc(1);
        });
    progress_bar.finish();
    let mut res: Vec<Output> = receiver.iter().collect();
    // parallel iteration messes up the order of scaffold ID's, so we fix that here.
    // however, within ID's windows should be ordered.
    res.sort_by_key(|x| x.id.clone());

    // write files
    // I pass over &res four times here...
    // pretty inefficient.
    // I should really separate these functions from main. It's just lazy.

    eprintln!("[+]\tWriting output to files");
    for i in &res {
        writeln!(
            window_file,
            // :.3 three decimal places
            "{}\t{}\t{}\t{}\t{:.3}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{:.3}\t{:.3}",
            i.id,
            i.start,
            i.end,
            i.gc_proportion,
            i.gc_skew,
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
    window_file.flush().unwrap();

    // these are the arrays, tab separated bed-like format.
    // TODO: add headers for each of these
    match kmer_maps.as_slice() {
        [two, three, four] => {
            // headers for dinucs
            let mut dinuc_headers = Vec::new();
            for key in two.map.keys().sorted() {
                dinuc_headers.push(key)
            }
            writeln!(
                window_file_2,
                "ID\tstart\tend\t{}",
                kmer_maps::WriteKmerValues(dinuc_headers)
            )
            .unwrap_or_else(|_| eprintln!("[-]\tError in writing to file."));

            for i in &res {
                writeln!(
                    window_file_2,
                    "{}\t{}\t{}\t{}",
                    i.id,
                    i.start,
                    i.end,
                    kmer_maps::WriteArray(i.divalues.clone())
                )
                .unwrap_or_else(|_| eprintln!("[-]\tError in writing to file."));
            }
            window_file_2.flush().unwrap();

            // headers for trinucs
            let mut trinuc_headers = Vec::new();
            for key in three.map.keys().sorted() {
                trinuc_headers.push(key)
            }
            writeln!(
                window_file_3,
                "ID\tstart\tend\t{}",
                kmer_maps::WriteKmerValues(trinuc_headers)
            )
            .unwrap_or_else(|_| eprintln!("[-]\tError in writing to file."));
            for i in &res {
                writeln!(
                    window_file_3,
                    "{}\t{}\t{}\t{}",
                    i.id,
                    i.start,
                    i.end,
                    kmer_maps::WriteArray(i.trivalues.clone())
                )
                .unwrap_or_else(|_| eprintln!("[-]\tError in writing to file."));
            }
            window_file_3.flush().unwrap();
            //
            // headers for tetranucs
            let mut tetranuc_headers = Vec::new();
            for key in four.map.keys().sorted() {
                tetranuc_headers.push(key)
            }
            writeln!(
                window_file_4,
                "ID\tstart\tend\t{}",
                kmer_maps::WriteKmerValues(tetranuc_headers)
            )
            .unwrap_or_else(|_| eprintln!("[-]\tError in writing to file."));
            for i in &res {
                writeln!(
                    window_file_4,
                    "{}\t{}\t{}\t{}",
                    i.id,
                    i.start,
                    i.end,
                    kmer_maps::WriteArray(i.tetravalues.clone())
                )
                .unwrap_or_else(|_| eprintln!("[-]\tError in writing to file."));
            }
            window_file_4.flush().unwrap();
        }
        [..] => {} // Needed to make the patterns exhaustive
    }

    eprintln!("[+]\tOutput written to directory: ./fw_out/{}", output);
}
