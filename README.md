# fasta_windows

Fast statistics in windows over a genome in fasta format.
- GC content
- GC proportion
- GC skew
- Shannon entropy
- Proportion of G's, C's, A's, T's, N's
- Di/tri/tetranucleotide shannon diversity

## Usage

Fewer options than previous versions, as di/tri/tetranucleotide diversity is calculated instead of user input for kmer length.

```
Fasta windows 0.2.0
Max Brown <mb39@sanger.ac.uk>
Quickly compute statistics over a fasta file in windows.

USAGE:
    fasta_windows [OPTIONS] --fasta <fasta> --output <output>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -c, --canonical_kmers <canonical_kmers>    Should the canonical kmers be calculated? Boolean, input true or false.
                                               [default: false]
    -f, --fasta <fasta>                        The input fasta file.
    -o, --output <output>                      Output filename for the CSV (without extension).
    -w, --window_size <window_size>            Integer size of window for statistics to be computed over. [default:
                                               1000]
```

You have to complile yourself. <a href="https://www.rust-lang.org/tools/install">Download rust</a>, clone this repo, and then run:

`cargo build --release`

Compiling may take a couple of minutes. This will then make the compiled binary in the `target/release` directory.

Run `./target/release/fasta_windows --help` to display the help message in the terminal.

For example, to iterate over a fasta file in windows of 100 base pairs, computing trinucleotide diversity:

`./target/release/fasta_windows --fasta /path/to/your/fasta --kmer_size 3 --window_size 100 --output /path/to/output`

Or to use default kmer length and windows (1kb), and calculate kmer count distance from window to genome wide:

`./target/release/fasta_windows --fasta /path/to/your/fasta --output /path/to/output --kmer_distance true`

## Output & benchmarks

### Output

Output is now a tsv with bed-like format in the first three columns:



### Updates & bugs

Please use, test, and let me know if there are any bugs or features you want implemented. Either raise an issue, or email me (see email in usage).