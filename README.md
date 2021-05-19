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

The default window size is 1kb.

## Output & benchmarks

The only annoying overhead at the moment is

### Output

Output is now a tsv with bed-like format in the first three columns:

```
ID      start   end     GC_prop GC_skew Shannon_entropy Prop_Gs Prop_Cs Prop_As Prop_Ts Prop_Ns Dinucleotide_Shannon_false      Trinucleotide_Shannon_false Tetranucleotide_Shannon_false
SUPER_1 0       1000    0.452   -0.2699115      1.928739902650348       0.165   0.287   0.361   0.187   0       2.6459008551823886 3.928519261697192        5.133591371839395
SUPER_1 1000    2000    0.34    -0.33529413     1.8955852733798557      0.113   0.227   0.346   0.314   0       2.6167836230348853 3.8722623719711  5.015274395434933
SUPER_1 2000    3000    0.388   -0.91237116     1.627180642639534       0.017   0.371   0.407   0.205   0       1.858410057857901  2.0494842744481336       2.09550303360082
SUPER_1 3000    4000    0.634   -0.16719243     1.9326861804290671      0.264   0.37    0.199   0.167   0       2.6709080342328937 3.9796052529877928       5.214642263562323
SUPER_1 4000    5000    0.591   -0.18443316     1.9543596224588031      0.241   0.35    0.236   0.173   0       2.701288242079077  4.0199349099815  5.232032920693032
SUPER_1 5000    6000    0.599   -0.22871453     1.9477765017208162      0.231   0.368   0.212   0.189   0       2.6791744546822502 3.990975528462955        5.20873424760944
SUPER_1 6000    7000    0.596   -0.16442953     1.9605365300597528      0.249   0.347   0.214   0.19    0       2.6935889794270693 3.9940001045093587       5.206001722737892
SUPER_1 7000    8000    0.602   -0.19269103     1.9503405864559629      0.243   0.359   0.178   0.22    0       2.671998818221988  3.9740681661842774       5.184128931560171
SUPER_1 8000    9000    0.453   -0.21412803     1.9767106890447885      0.178   0.275   0.292   0.255   0       2.7253730593872803 4.030994655335826        5.2367638611178435
```


### Updates & bugs

Please use, test, and let me know if there are any bugs or features you want implemented. Either raise an issue, or email me (see email in usage).