# fasta_windows.rs

Fast statistics in windows over a genome in fasta format.
- GC content
- GC skew
- Kmer diversity (canonical or not; tested kmers up to 31 with little performance dip)
- Kmer distance (optional), the euclidean distance of the kmer profile of the current window to the reference.
- Shannon entropy

## Usage

```
Fasta windows 0.1.2
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
    -d, --kmer_distance <kmer_distance>        Calculate kmer count distance to reference? Boolean, input true or false.
                                               NOTE experimental, needs QC. [default: false]
    -k, --kmer_size <kmer_size>                Size of kmer to determine the diversity of in windows. [default: 4]
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

Or to use default kmer length and windows, and calculate kmer count distance from window to genome wide:

`./target/release/fasta_windows --fasta /path/to/your/fasta --output /path/to/output --kmer_distance true`

## Output & benchmarks

### Output

Output is a CSV file with headers:

```
ID,window,GC_percent,GC_skew,Shannon_entropy,4mer_diversity_canonical_false
NC_003070.9,1000,32.7,-0.10703364,1.8933459501670153,209
NC_003070.9,2000,31.7,-0.05362776,1.8994159317068182,206
NC_003070.9,3000,32.5,-0.1323077,1.905439536788656,216
NC_003070.9,4000,33.399998,-0.011976048,1.9184021842094099,231
NC_003070.9,5000,39.5,0.07848101,1.965667840754946,236
NC_003070.9,6000,36.199997,0.13812155,1.9389668067700883,229
NC_003070.9,7000,38,0.02631579,1.9565065330018037,233
NC_003070.9,8000,34.1,-0.06158358,1.9070730536871001,217
NC_003070.9,9000,34.4,-0.046511628,1.9220237551539758,227
```

And also currently printed to stdout is the number of sequences in the fasta file (ideally chromosomes), length of the total genome, and the N50.

```
[+]	NC_003070.9 processed.
[+]	NC_003071.7 processed.
[+]	NC_000932.1 processed.
[+]	NC_003074.8 processed.
[+]	NC_003075.7 processed.
[+]	NC_003076.8 processed.
[+]	NC_037304.1 processed.
[+]	Global stats:
Number of contigs/chromosomes: 7
    Total length of genome: 119668634
    The N50 of this genome: 23459830
```

### Tests 

Operating on the concatenated genome of *Arabidopsis thaliana* including plastid & mitochondrial genomes (~120Mb). Source fastas <a href="https://www.ncbi.nlm.nih.gov/genome/?term=arabidopsis%20thaliana">here</a>. Calculating canonical kmers takes slightly longer, as does kmer count distance to reference. The latter requires two passes of the genome. Testing done on my Mac, so perhaps a pinch of salt.

Example command:

`time ./target/release/fasta_windows --fasta Athaliana_1_5_m_c.fasta --output test`

```
real	0m17.823s
user    0m17.068s
sys 0m0.614s
```
