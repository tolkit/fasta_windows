// external imports
// std
use std::fs::{File, create_dir_all};
use std::io::LineWriter;
use std::io::prelude::*;
use std::collections::HashMap;
// non-std
extern crate clap; // forgot why I needed extern crate.
use clap::{App, Arg, value_t};
use bio::io::fasta;
// internal imports
use fasta_windows::gcu8::gcu8;
use fasta_windows::kmeru8::kmeru8;
use fasta_windows::wgs::wgs;
use fasta_windows::utils::utils;

// TODO: can I implement multiple threads?
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
                 .help("Should the canonical kmers be calculated? Boolean, input true or false.")
                 .takes_value(true)
                 .default_value("false"))
        .arg(Arg::with_name("kmer_distance")
                 .short("d")
                 .long("kmer_distance")
                 .help("Calculate kmer count distance to reference? Boolean, input true or false. NOTE experimental, needs QC.")
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
    let kmer_distance = value_t!(matches.value_of("kmer_distance"), bool).unwrap_or_else(|e| e.exit());

    // create directory for output
    create_dir_all("./fw_out/");

    // initiate the output CSV
    let output_file = format!("./fw_out/{}{}", output, ".csv");
    let file = File::create(&output_file).unwrap();
    let mut file = LineWriter::new(file);
    // and write the headers
    if kmer_distance {
        writeln!(file, "ID,window,GC_percent,GC_skew,{}mer_diversity_canonical_{},kmer_distance", kmer_size, canonical_kmers).unwrap();
    } else {
        writeln!(file, "ID,window,GC_percent,GC_skew,{}mer_diversity_canonical_{}", kmer_size, canonical_kmers).unwrap();
    }
        
    // read in the fasta from file
    let reader = fasta::Reader::from_file(input_fasta).expect("[-]\tPath invalid.");

    // the sequence lengths of each fasta record.
    let mut read_lengths: Vec<usize> = Vec::new();
    // get kmer hash of entire genome
    let kmer_hash = &mut HashMap::new();

    if kmer_distance {
        // creating new fasta reader here shouldnt have large overhead.
        let kmer_reader = fasta::Reader::from_file(input_fasta).expect("[-] Path invalid.");
        println!("[+]\tFirst pass of genome.");
        for result in kmer_reader.records() {
            let record = result.expect("Error during fasta record parsing");
            // the current kmer hash 
            let kmer_current = kmeru8::kmer_diversity(record.seq(), kmer_size, canonical_kmers).kmer_hash;
            // merge current kmer hash with previous iteration
            let kmer_hash = utils::merge_hashmap_ip(kmer_hash, kmer_current);
        }
        println!("[+]\tWhole genome kmer HashMap made.")
    }

    // iterate over fasta records
    if kmer_distance {
        println!("[+]\tSecond pass of genome.");
    }
    for result in reader.records() {
        let record = result.expect("[-]\tError during fasta record parsing.");

        // for the stats at the end.
        read_lengths.push(record.seq().len());
        // initiate a counter for the windows
        let mut counter = window_size;
        // begin sliding windows
        let windows = record.seq().chunks(window_size);

        for win in windows {
            let win_gc = gcu8::gc_content(win);
            let kmer_stats = kmeru8::kmer_diversity(win, kmer_size, canonical_kmers);
            let kmer_distance_value = utils::create_kmer_distance(kmer_stats.kmer_hash, kmer_hash);
            // ugly way of handling this but...
            if kmer_distance {
                writeln!(file, "{},{},{},{},{},{}", record.id(), counter, win_gc.gc_content, win_gc.gc_skew, kmer_stats.kmer_diversity, kmer_distance_value).unwrap();    
            } else {
                writeln!(file, "{},{},{},{},{}", record.id(), counter, win_gc.gc_content, win_gc.gc_skew, kmer_stats.kmer_diversity).unwrap();
            }
            // re-set the counter if counter > length of current sequence
            if counter < record.seq().len() {
                counter += window_size
            } else {
                counter = 0
            }
        }
        println!("[+]\t{} processed.", record.id());
    }

    println!("[+]\tGlobal stats:");
    // eventually will write to file.
    let genome_stats = wgs::collect_genome_stats(read_lengths);
    println!("\tNumber of contigs/chromosomes: {}", genome_stats.no_reads);
    println!("\tTotal length of genome: {}", genome_stats.genome_length);
    println!("\tThe N50 of this genome: {}", genome_stats.n50);
}