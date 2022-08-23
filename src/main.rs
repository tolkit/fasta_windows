// Max Brown & Matthieu Muffato 2022; Wellcome Sanger Institute

// std imports
use std::fs::{create_dir_all, File};
use std::io::BufWriter;

// non-std
use anyhow::{Context, Result};
use clap::{Arg, Command};
// internal imports
use fasta_windows::fw::fasta_windows;

fn main() -> Result<()> {
    // command line options
    let matches = Command::new("Fasta windows")
        .version(clap::crate_version!())
        .arg_required_else_help(true)
        .author("Max Brown <mb39@sanger.ac.uk>")
        .about("Quickly compute statistics over a fasta file in windows.")
        .arg(
            Arg::new("fasta")
                .short('f')
                .long("fasta")
                .takes_value(true)
                .required(true)
                .help("The input fasta file."),
        )
        .arg(
            Arg::new("window_size")
                .short('w')
                .long("window_size")
                .help("Integer size of window for statistics to be computed over.")
                .takes_value(true)
                .default_value("1000"),
        )
        .arg(
            Arg::new("description")
                .short('d')
                .long("description")
                .help("Add an extra column to _windows.tsv output with fasta header descriptions."),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .help("Output filename for the TSV's (without extension).")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::new("masked")
                .short('m')
                .long("masked")
                .help("Consider only uppercase nucleotides in the calculations."),
        )
        .get_matches();
    // parse command line options
    let output = matches
        .value_of("output")
        .context("Could not find output in CLI")?;

    // create directory for output
    if let Err(e) = create_dir_all("./fw_out/") {
        eprintln!("[-]\tCreate directory error: {}", e.to_string());
    }

    // initiate the output TSV files
    let output_file_0 = format!("./fw_out/{}{}", output, "_freq_windows.tsv");
    let window_file_0 = File::create(&output_file_0)?;
    let window_file_0 = BufWriter::new(window_file_0);

    let output_file_1 = format!("./fw_out/{}{}", output, "_mononuc_windows.tsv");
    let window_file_1 = File::create(&output_file_1)?;
    let window_file_1 = BufWriter::new(window_file_1);

    let output_file_2 = format!("./fw_out/{}{}", output, "_dinuc_windows.tsv");
    let window_file_2 = File::create(&output_file_2)?;
    let window_file_2 = BufWriter::new(window_file_2);

    let output_file_3 = format!("./fw_out/{}{}", output, "_trinuc_windows.tsv");
    let window_file_3 = File::create(&output_file_3)?;
    let window_file_3 = BufWriter::new(window_file_3);

    let output_file_4 = format!("./fw_out/{}{}", output, "_tetranuc_windows.tsv");
    let window_file_4 = File::create(&output_file_4)?;
    let window_file_4 = BufWriter::new(window_file_4);

    // pass this function the matches from clap && the four files to write.
    fasta_windows(
        &matches,
        window_file_0,
        window_file_1,
        window_file_2,
        window_file_3,
        window_file_4,
    )?;

    Ok(())
}
