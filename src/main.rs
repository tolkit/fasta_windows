// Max Brown, Matthieu Muffato, and Rich Challis 2023
// Wellcome Sanger Institute

use std::fs::{create_dir_all, File};
use std::io::BufWriter;
use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::{crate_authors, value_parser, Arg, ArgAction, Command};
use fasta_windows::fw::fasta_windows;

fn main() -> Result<()> {
    let cmd = Command::new("Fasta windows")
        .version(clap::crate_version!())
        .arg_required_else_help(true)
        .author(crate_authors!())
        .about("Quickly compute statistics over a fasta file in windows.")
        .arg(
            Arg::new("fasta")
                .short('f')
                .long("fasta")
                .num_args(1)
                .required(true)
                .value_parser(value_parser!(PathBuf))
                .help("The input fasta file."),
        )
        .arg(
            Arg::new("window_size")
                .short('w')
                .long("window_size")
                .help("Integer size of window for statistics to be computed over.")
                .num_args(1)
                .value_parser(value_parser!(usize))
                .default_value("1000"),
        )
        .arg(
            Arg::new("description")
                .short('d')
                .long("description")
                .action(ArgAction::SetTrue)
                .help("Add an extra column to _windows.tsv output with fasta header descriptions."),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .help("Output filename for the TSV's (without extension).")
                .value_parser(value_parser!(PathBuf))
                .num_args(1)
                .required(true),
        )
        .arg(
            Arg::new("masked")
                .short('m')
                .long("masked")
                .action(ArgAction::SetTrue)
                .help("Consider only uppercase nucleotides in the calculations."),
        )
        .arg(
            Arg::new("ctw")
                .short('c')
                .long("ctw")
                .action(ArgAction::SetTrue)
                .help("Calculate the Context-Tree Weighting (how compressible a sequence is)."),
        );

    #[cfg(feature = "entropy")]
    let cmd = cmd.arg(
        Arg::new("entropy")
            .short('e')
            .long("entropy")
            .action(ArgAction::SetTrue)
            .help(
                "Entropy mode: output a single BED file of Shannon entropy per window. \
                 Skips all k-mer and nucleotide-composition computation for maximum speed.",
            ),
    );

    let matches = cmd.get_matches();

    let output = matches
        .get_one::<PathBuf>("output")
        .context("Could not find output in CLI")?
        .display();

    if let Err(e) = create_dir_all("./fw_out/") {
        eprintln!("[-]\tCreate directory error: {}", e);
    }

    #[cfg(feature = "entropy")]
    if matches.get_flag("entropy") {
        use fasta_windows::entropy::entropy_windows;
        let bed_file = BufWriter::new(File::create(format!("./fw_out/{output}_entropy.bed"))?);
        return entropy_windows(&matches, bed_file);
    }

    let window_file_0 =
        BufWriter::new(File::create(format!("./fw_out/{output}_freq_windows.tsv"))?);
    let window_file_1 = BufWriter::new(File::create(format!(
        "./fw_out/{output}_mononuc_windows.tsv"
    ))?);
    let window_file_2 = BufWriter::new(File::create(format!(
        "./fw_out/{output}_dinuc_windows.tsv"
    ))?);
    let window_file_3 = BufWriter::new(File::create(format!(
        "./fw_out/{output}_trinuc_windows.tsv"
    ))?);
    let window_file_4 = BufWriter::new(File::create(format!(
        "./fw_out/{output}_tetranuc_windows.tsv"
    ))?);

    fasta_windows(
        &matches,
        window_file_0,
        window_file_1,
        window_file_2,
        window_file_3,
        window_file_4,
    )
}
