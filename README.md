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

```
ID  start   end GC_content  GC_prop GC_skew Shannon_entropy Prop_Gs Prop_Cs Prop_As Prop_Ts Prop_Ns Dinucleotide_Shannon_false  Trinucleotide_Shannon_false Tetranucleotide_Shannon_false
   NC_000932.1 0   1000    35.100002   0.351   -0.08262108 1.9262736915905785  0.161   0.19    0.364   0.285   0   2.657852221603254   3.959546193744067   5.178306365245311
   NC_000932.1 1000    2000    36  0.36    -0.03888889 1.933130471343546   0.173   0.187   0.365   0.275   0   2.6705953385790915  3.9681692272585343  5.1957269237412795
   NC_000932.1 2000    3000    31.5    0.315   -0.07936508 1.8951002428635593  0.145   0.17    0.366   0.319   0   2.6045409546215144  3.8710188243791652  5.047995206208306
   NC_000932.1 3000    4000    30.4    0.304   -0.019736841    1.8842108651215588  0.149   0.155   0.369   0.327   0   2.6000265084811076  3.8757372067753177  5.064073968458029
   NC_000932.1 4000    5000    26.199999   0.262   -0.015267176    1.829701464570813   0.129   0.133   0.37    0.368   0   2.52333440322612    3.7462149281116606  4.871866243052142
    NC_000932.1 5000    6000    31.3    0.313   -0.15654951 1.8892838217462478  0.132   0.181   0.323   0.364   0   2.593172762071583   3.853857537428029   5.028882758920603
   NC_000932.1 6000    7000    28.5    0.285   -0.024561403    1.861928417304251   0.139   0.146   0.363   0.352   0   2.5714885317604574  3.8324913747638694  4.996857758056243
   NC_000932.1 7000    8000    29.9    0.299   -0.12374582 1.8712601936539917  0.131   0.168   0.314   0.387   0   2.5692697851487956  3.8126901042473804  4.973042839086266
  NC_000932.1 8000    9000    23.199999   0.232   -0.03448276 1.7800647412527244  0.112   0.12    0.366   0.402   0   2.446987915911496   3.6225846868349842  4.689852544893519
  NC_000932.1 9000    10000   29.4    0.294   -0.020408163    1.8735437994390223  0.144   0.15    0.36    0.346   0   2.5735874748973933  3.818744232935321   4.973867185531613
  NC_000932.1 10000   11000   40.6    0.406   0.02955665  1.971743246519977   0.209   0.197   0.275   0.319   0   2.7263411622643314  4.039131168815778   5.259445686803222
```


### Updates & bugs

Please use, test, and let me know if there are any bugs or features you want implemented. Either raise an issue, or email me (see email in usage).