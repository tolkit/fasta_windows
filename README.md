# fasta_windows.rs

A re-write of the fasta_window_stats written in Rust.

Currently only GC content and kmer diversity in sliding windows are implemented.

## Usage

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
    -c, --canonical_kmers <canonical_kmers>    Should the canonical kmers be calculated? Bool, input true or false. [default: false]
    -f, --fasta <fasta>                        The input fasta file.
    -k, --kmer_size <kmer_size>                Size of kmer to determine the diversity of in windows. [default: 4]
    -o, --output <output>                      Output filename for the CSV (without extension).
    -w, --window_size <window_size>            Integer size of window for statistics to be computed over. [default: 1000]
```

I think you have to complile yourself, the binaries I've made are quite big. You will need to <a href="https://www.rust-lang.org/tools/install">download rust</a>, clone this repo, and then run:

`cargo build --release`

This will then make the compiled binary in the `target/release` directory.

Run `./target/release/fasta_windows --help` to display the help message in the terminal.

For example, to iterate over a fasta file in windows of 100 base pairs, computing trinucleotide diversity:

`./target/release/fasta_windows --fasta /path/to/your/fasta --kmer_size 3 --window_size 100 --output /path/to/output`

## Output & benchmarks

Output is a CSV file with headers:

```
ID,window,GC_percent,GC_skew,4mer_diversity_canonical_true
MK070895.1,1000,36.199997,-0.049723756,132
MK070895.1,2000,37.5,-0.056,135
MK070895.1,3000,34.3,-0.0670554,129
MK070895.1,4000,32.8,-0.030487806,129
MK070895.1,5000,31.5,0.034920637,130
MK070895.1,6000,34.4,-0.034883723,125
MK070895.1,7000,27.7,-0.010830325,128
MK070895.1,8000,31.600002,-0.07594936,132
MK070895.1,9000,30.199999,-0.13245033,129
```

### Tests 

Operating on the concatenated genome of *Arabidopsis thaliana* including plastid & mitochondrial genomes... Source fastas <a href="https://www.ncbi.nlm.nih.gov/genome/?term=arabidopsis%20thaliana">here</a>.

Command:

`time ./target/release/fasta_windows --fasta Athaliana_1_5_m_c.fasta --output test`

```
NC_003070.9 processed.
NC_003071.7 processed.
NC_000932.1 processed.
NC_003074.8 processed.
NC_003075.7 processed.
NC_003076.8 processed.
NC_037304.1 processed.
--------------------------
Number of contigs/chromosomes: 7
Total length of genome: 119668634
The N50 of this genome: 23459830

real	0m19.239s
user	0m18.516s
sys	0m0.622s
```

Lost a little time in computing canonical kmers, looking to implement in u8...