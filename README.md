# fasta_windows

<p align="center">
    <img width="300" height="132" src="https://www.darwintreeoflife.org/wp-content/themes/dtol/dist/assets/gfx/dtol-logo-round.png">
</p>

Written for Darwin Tree of Life chromosomal level genome assemblies. The executable takes a fasta formatted file and calculates some statistics of interest in windows:

- GC content
- GC proportion
- GC skew
- Proportion of G's, C's, A's, T's, N's
- Shannon entropy
- Di/tri/tetranucleotide shannon diversity
- Di/tri/tetranucleotide frequency arrays

Output files can be visualised using <a href="https://github.com/tolkit/fw_plot">fw_plot</a> or grouped using <a href="https://github.com/tolkit/fw_group">fw_group</a>.

## Download

The easiest way to get `fasta_windows` is through <b><a href="https://anaconda.org/bioconda/fasta_windows">conda/bioconda</a></b>.

```bash
conda create -n fasta_windows -c bioconda fasta_windows
```

## Usage

```
Fasta windows 0.2.3
Max Brown <mb39@sanger.ac.uk>
Quickly compute statistics over a fasta file in windows.

USAGE:
    fasta_windows [FLAGS] [OPTIONS] --fasta <fasta> --output <output>

FLAGS:
    -d, --description    Add an extra column to _windows.tsv output with fasta header descriptions.
    -h, --help           Prints help information
    -m, --masked         Consider only uppercase nucleotides in the calculations.
    -V, --version        Prints version information

OPTIONS:
    -f, --fasta <fasta>                The input fasta file.
    -o, --output <output>              Output filename for the TSV's (without extension).
    -w, --window_size <window_size>    Integer size of window for statistics to be computed over. [default: 1000]
```

## Building

Building <a href="https://www.rust-lang.org/tools/install">requires Rust</a>. 

```bash
git clone https://github.com/tolkit/fasta_windows
cd fasta_windows
cargo build --release
# ./target/release/fasta_windows is the executable
# show help
./target/release/fasta_windows --help
```

The default window size is 1kb.

## Output

Output is now a tsv with bed-like format in the first three columns:

```
ID      start   end     GC_prop GC_skew Shannon_entropy Prop_Gs Prop_Cs Prop_As Prop_Ts Prop_Ns Dinucleotide_Shannon_false      Trinucleotide_Shannon_false Tetranucleotide_Shannon_false
SUPER_1 0       1000    0.452   -0.270  1.929   0.165   0.287   0.361   0.187   0       2.646   3.929   5.134
SUPER_1 1000    2000    0.34    -0.335  1.896   0.113   0.227   0.346   0.314   0       2.617   3.872   5.015
SUPER_1 2000    3000    0.388   -0.912  1.627   0.017   0.371   0.407   0.205   0       1.858   2.049   2.096
SUPER_1 3000    4000    0.634   -0.167  1.933   0.264   0.37    0.199   0.167   0       2.671   3.980   5.215
SUPER_1 4000    5000    0.591   -0.184  1.954   0.241   0.35    0.236   0.173   0       2.701   4.020   5.232
SUPER_1 5000    6000    0.599   -0.229  1.948   0.231   0.368   0.212   0.189   0       2.679   3.991   5.209
SUPER_1 6000    7000    0.596   -0.164  1.961   0.249   0.347   0.214   0.19    0       2.694   3.994   5.206
SUPER_1 7000    8000    0.602   -0.193  1.950   0.243   0.359   0.178   0.22    0       2.672   3.974   5.184
SUPER_1 8000    9000    0.453   -0.214  1.977   0.178   0.275   0.292   0.255   0       2.725   4.031   5.237
```

Also output (non-optional at the moment), are three more TSV's, which are the arrays of di/tri/tetranucleotide frequencies in each window. These files are large, especially as tetranucleotide frequencies will contain 4e4 columns. The kmers are sorted lexicographically from left -> right (AA(AA) to TT(TT)).

e.g. for dinucleotide frequencies:

```
ID	start	end	AA	AC	AG	AT	CA	CC	CG	CT	GA	GC	GG	GT	TA	TC	TG	TT
SUPER_1 0       1000    122     120     45      73      134     68      39      46      50      55      45      15      54      44 36       53
SUPER_1 1000    2000    140     83      32      90      85      54      22      66      30      25      19      39      91      65 40       118
SUPER_1 2000    3000    216     181     4       5       4       181     5       181     3       8       3       3       183     1  516
SUPER_1 3000    4000    40      61      54      44      80      137     86      66      54      99      76      35      24      73 48       22
SUPER_1 4000    5000    55      68      75      38      88      138     66      57      58      78      59      46      35      65 41       32
SUPER_1 5000    6000    32      71      63      46      85      137     71      75      65      66      65      34      30      94 31       34
SUPER_1 6000    7000    47      62      63      42      91      132     60      64      58      84      74      32      18      69 51       52
SUPER_1 7000    8000    29      49      64      35      67      143     52      97      58      82      72      31      24      85 55       56
SUPER_1 8000    9000    114     67      43      68      63      86      52      73      51      49      43      35      64      73 40       78
SUPER_1 9000    10000   97      97      44      63      72      95      50      67      46      44      33      46      85      49 42       69
```

### Comments, updates & bugs

As of version 0.2.2, I've removed canonical kmers as an option; it was really computationally expensive and I couldn't think of a way to efficienty add it in. End users that wish this are pointed in the direction of <a href="https://github.com/tolkit/fw_group">fw_group</a>, which will at some point soon provide this functionality.

The masked (-m) flag only affects GC content, GC proportion, GC skew, proportion of G's, C's, A's, T's, N's. Kmers are coerced to uppercase automatically. Shannon index counts only uppercase nucleotides.

Please use, test, and let me know if there are any bugs or features you want implemented. Either raise an issue, or email me (see email in usage).