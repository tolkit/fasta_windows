// external imports
use std::str;
use std::fs::File;
use std::io::LineWriter;
use std::io::prelude::*;
extern crate clap; // forgot why I needed extern crate.
use clap::{App, Arg, value_t};
use bio::io::fasta;

// internal imports
use fasta_windows::gc::gc;
use fasta_windows::windows::windows;
use fasta_windows::kmer::kmer;

// TODO: can I implement multiple threads?

fn main() {
    // command line options
    let matches = App::new("Fasta windows")
        .version("0.1.0")
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

    // initiate the output CSV
    let output_file = format!("{}{}", output, ".csv");
    let file = File::create(&output_file).unwrap();
    let mut file = LineWriter::new(file);
    // and write the headers
    writeln!(file, "ID,window,GC_percent,kmer_diversity").unwrap();
        
    // read in the fasta from file
    let reader = fasta::Reader::from_file(input_fasta).expect("Path invalid.");

    let mut nb_reads = 0;
    let mut nb_bases = 0;

    // iterate over fasta records
    for result in reader.records() {
        let record = result.expect("Error during fasta record parsing");

        nb_reads += 1;
        nb_bases += record.seq().len();

        // https://stackoverflow.com/questions/19076719/how-do-i-convert-a-vector-of-bytes-u8-to-a-string
        // the process of turning u8 into str appears to have very little overhead
        // and is more human readable? TODO test this.

        let nucleotide_string = match str::from_utf8(record.seq()) {
            Ok(v) => v,
            Err(e) => panic!("Invalid UTF-8 sequence: {}", e),
        };

        let mut counter = window_size;

        // here sliding windows
        let windows = windows::char_windows(nucleotide_string, window_size, window_size);

        for win in windows {
            let win_gc = gc::gc_content(win);
            let no_kmers = kmer::kmer_diversity(win, kmer_size);

            // not sure what the overhead of this is compared to writing at the end
            writeln!(file, "{},{},{},{}", record.id(), counter, win_gc, no_kmers).unwrap();
            
            // re-set the counter if 
            if counter < record.seq().len() {
                counter += window_size
            } else {
                counter = 0
            }
        }

        let gc = gc::gc_content(nucleotide_string);
        println!("{} processed.", record.id());
        println!("Total GC: {}", gc);

    }

    println!("Number of contigs/chromosomes: {}", nb_reads);
    println!("Number of bases processed: {}", nb_bases);
}