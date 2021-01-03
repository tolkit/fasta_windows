# fasta_windows.rs

A re-write of the fasta_window_stats written in Rust.

Currently only GC content and kmer diversity in sliding windows are implemented.

## Usage

The compiled binary is present in `./target/release/fasta_windows`. Usage is as follows:

```
Fasta windows 0.1.0
Max Brown <mb39@sanger.ac.uk>
Quickly compute statistics over a fasta file in windows.

USAGE:
    fasta_windows [OPTIONS] --fasta <fasta> --output <output>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -f, --fasta <fasta>                The input fasta file.
    -k, --kmer_size <kmer_size>        Size of kmer to determine the diversity of in windows. [default: 4]
    -o, --output <output>              Output filename for the CSV (without extension).
    -w, --window_size <window_size>    Integer size of window for statistics to be computed over. [default: 1000]
```

Run `./target/release/fasta_windows --help` to display this message in the terminal.

For example, to iterate over a fasta file in windows of 100 base pairs, computing trinucleotide diversity:

`./target/release/fasta_windows --fasta /path/to/your/fasta --kmer_size 3 --window_size 100 --output /path/to/output`

## Output & benchmarks

Output is a CSV file with headers:

```
ID,window,GC_percent,kmer_diversity
...
...
```

Operating on the concatenated genome of *Arabidopsis thaliana* including plastid & mitochondrial genomes, takes only 8 seconds. Source fastas <a href="https://www.ncbi.nlm.nih.gov/genome/?term=arabidopsis%20thaliana">here</a>.

Printed to stout (using time prefix):

```
NC_003070.9 processed.
Total GC: 0.35679775
NC_003071.7 processed.
Total GC: 0.3585966
NC_000932.1 processed.
Total GC: 0.3629384
NC_003074.8 processed.
Total GC: 0.3632182
NC_003075.7 processed.
Total GC: 0.36198115
NC_003076.8 processed.
Total GC: 0.35925233
NC_037304.1 processed.
Total GC: 0.4479212
Number of contigs/chromosomes: 7
Number of bases processed: 119668634

real	0m8.845s
user	0m8.468s
sys	0m0.365s
```