# fasta_windows.rs

A re-write of the fasta_window_stats written in Rust.

Currently only GC content and kmer diversity in sliding windows are implemented.

## Usage

```
Fasta windows 0.1.1
Max Brown <mb39@sanger.ac.uk>
Quickly compute statistics over a fasta file in windows.

USAGE:
    fasta_windows [OPTIONS] --fasta <fasta> --output <output>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -c, --canonical_kmers <canonical_kmers>    Should the canonical kmers be calculated? Bool, input true or false.
                                               [default: false]
    -f, --fasta <fasta>                        The input fasta file.
    -k, --kmer_size <kmer_size>                Size of kmer to determine the diversity of in windows. [default: 4]
    -o, --output <output>                      Output filename for the CSV (without extension).
    -w, --window_size <window_size>            Integer size of window for statistics to be computed over. [default:
                                               1000]
```

You have to complile yourself. <a href="https://www.rust-lang.org/tools/install">Download rust</a>, clone this repo, and then run:

`cargo build --release`

This will then make the compiled binary in the `target/release` directory.

Run `./target/release/fasta_windows --help` to display the help message in the terminal.

For example, to iterate over a fasta file in windows of 100 base pairs, computing trinucleotide diversity:

`./target/release/fasta_windows --fasta /path/to/your/fasta --kmer_size 3 --window_size 100 --output /path/to/output`

## Output & benchmarks

### Output

Output is a CSV file with headers:

```
ID,window,GC_percent,GC_skew,4mer_diversity_canonical_false
NC_003070.9,1000,32.7,-0.10703364,209
NC_003070.9,2000,31.7,-0.05362776,206
NC_003070.9,3000,32.5,-0.1323077,216
NC_003070.9,4000,33.399998,-0.011976048,231
NC_003070.9,5000,39.5,0.07848101,236
NC_003070.9,6000,36.199997,0.13812155,229
NC_003070.9,7000,38,0.02631579,233
NC_003070.9,8000,34.1,-0.06158358,217
NC_003070.9,9000,34.4,-0.046511628,227
```

And also currently printed to stdout is the number of sequences in the fasta file (ideally chromosomes), length of the total genome, and the N50.

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
```

### Tests 

Operating on the concatenated genome of *Arabidopsis thaliana* including plastid & mitochondrial genomes (~120Mb). Source fastas <a href="https://www.ncbi.nlm.nih.gov/genome/?term=arabidopsis%20thaliana">here</a>.

Command:

`time ./target/release/fasta_windows --fasta Athaliana_1_5_m_c.fasta --output test`

```
real	0m16.211s
user	0m15.670s
sys	0m0.480s
```
