// external imports
use std::fs::File;
use std::io::LineWriter;
use std::io::prelude::*;
extern crate clap; // forgot why I needed extern crate.
use clap::{App, Arg, value_t};
use bio::io::fasta;

// internal imports
use fasta_windows::gcu8::gcu8;
use fasta_windows::kmeru8::kmeru8;
use fasta_windows::wgs::wgs;

// TODO: can I implement multiple threads?
//     : compare kmer distribution per window to genome wide distribution
//     : separate file for genome stats.

fn main() {
    // command line options
    let matches = App::new("Fasta windows")
        .version(clap::crate_version!())
        .author("Max Brown <mb39@sanger.ac.uk>")
        .about("Quickly compute statistics over a fasta file in windows.")
        .arg(Arg::with_name("fasta")
                 .short("f")
                 .long("fasta")
                 .takes_value(true)
                 .required(true)
                 .help("The input fasta file."))
        .arg(Arg::with_name("window_size")
                 .short("w")
                 .long("window_size")
                 .help("Integer size of window for statistics to be computed over.")
                 .takes_value(true)
                 .default_value("1000"))
        .arg(Arg::with_name("kmer_size")
                 .short("k")
                 .long("kmer_size")
                 .help("Size of kmer to determine the diversity of in windows.")
                 .takes_value(true)
                 .default_value("4"))
        .arg(Arg::with_name("canonical_kmers")
                 .short("c")
                 .long("canonical_kmers")
                 .help("Should the canonical kmers be calculated? Bool, input true or false.")
                 .takes_value(true)
                 .default_value("false"))
        .arg(Arg::with_name("output")
                 .short("o")
                 .long("output")
                 .help("Output filename for the CSV (without extension).")
                 .takes_value(true)
                 .required(true))
        .get_matches();
    
    // parse command line options
    let input_fasta = matches.value_of("fasta").unwrap();
    let output = matches.value_of("output").unwrap();
    let window_size = value_t!(matches.value_of("window_size"), usize).unwrap_or_else(|e| e.exit());
    let kmer_size = value_t!(matches.value_of("kmer_size"), usize).unwrap_or_else(|e| e.exit());
    let canonical_kmers = value_t!(matches.value_of("canonical_kmers"), bool).unwrap_or_else(|e| e.exit());

    // initiate the output CSV
    let output_file = format!("{}{}", output, ".csv");
    let file = File::create(&output_file).unwrap();
    let mut file = LineWriter::new(file);
    // and write the headers
    writeln!(file, "ID,window,GC_percent,GC_skew,{}mer_diversity_canonical_{}", kmer_size, canonical_kmers).unwrap();
        
    // read in the fasta from file
    let reader = fasta::Reader::from_file(input_fasta).expect("Path invalid.");

    // the sequence lengths of each fasta record.
    let mut read_lengths: Vec<usize> = Vec::new();

    // when I get round to implementing comparing kmer dist to genome wide
    // will need two passes of the genome

    // iterate over fasta records
    for result in reader.records() {
        let record = result.expect("Error during fasta record parsing");
        // for the stats at the end.
        read_lengths.push(record.seq().len());
        // initiate a counter for the windows
        let mut counter = window_size;
        // begin sliding windows
        let windows = record.seq().chunks(window_size);

        for win in windows {
            let win_gc = gcu8::gc_content(win);
            let no_kmers = kmeru8::kmer_diversity(win, kmer_size, canonical_kmers);

            // not sure what the overhead of this is compared to writing at the end
            writeln!(file, "{},{},{},{},{}", record.id(), counter, win_gc.gc_content, win_gc.gc_skew, no_kmers.kmer_diversity).unwrap();

            // re-set the counter if counter > length of current sequence
            if counter < record.seq().len() {
                counter += window_size
            } else {
                counter = 0
            }
        }
        println!("{} processed.", record.id());
    }

    println!("--------------------------");
    // for now, just print these stats at the end.
    // eventually will write to file.
    let genome_stats = wgs::collect_genome_stats(read_lengths);
    println!("Number of contigs/chromosomes: {}", genome_stats.no_reads);
    println!("Total length of genome: {}", genome_stats.genome_length);
    println!("The N50 of this genome: {}", genome_stats.n50);
}