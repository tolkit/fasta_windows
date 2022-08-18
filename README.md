# fasta_windows

<p align="center">
    <img width="300" height="132" src="https://www.darwintreeoflife.org/wp-content/themes/dtol/dist/assets/gfx/dtol-logo-round.png">
</p>

Written for Darwin Tree of Life chromosomal level genome assemblies. The executable takes a fasta formatted file and calculates some statistics of interest in windows:

- GC content
- GC proportion
- GC and AT skew
- Proportion of G's, C's, A's, T's, N's, and CpG's
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
ID	start	end	GC_prop	GC_skew	AT_skew	Shannon_entropy	Prop_Gs	Prop_Cs	Prop_As	Prop_Ts	Prop_Ns	CpG_prop	Dinucleotide_Shannon	Trinucleotide_Shannon	Tetranucleotide_Shannon
OV656674.1	0	1000	0.533	-0.088	-0.006	1.994	0.243	0.290	0.232	0.235	0.000	0.070	3.963	5.825	7.474
OV656674.1	1000	2000	0.645	0.048	0.025	1.937	0.338	0.307	0.182	0.173	0.000	0.120	3.862	5.751	7.526
OV656674.1	2000	3000	0.579	0.022	0.012	1.982	0.296	0.283	0.213	0.208	0.000	0.106	3.940	5.871	7.653
OV656674.1	3000	4000	0.541	0.039	-0.020	1.994	0.281	0.260	0.225	0.234	0.000	0.081	3.980	5.926	7.763
OV656674.1	4000	5000	0.585	0.084	-0.075	1.974	0.317	0.268	0.192	0.223	0.000	0.104	3.917	5.801	7.568
OV656674.1	5000	6000	0.529	-0.096	-0.006	1.994	0.239	0.290	0.234	0.237	0.000	0.068	3.938	5.740	7.297
OV656674.1	6000	7000	0.576	-0.118	0.075	1.976	0.254	0.322	0.228	0.196	0.000	0.079	3.948	5.884	7.666
OV656674.1	7000	8000	0.526	-0.004	0.084	1.996	0.262	0.264	0.257	0.217	0.000	0.065	3.975	5.903	7.692
OV656674.1	8000	9000	0.430	-0.093	0.088	1.980	0.195	0.235	0.310	0.260	0.000	0.054	3.955	5.899	7.719
```

Also output (non-optional at the moment), are three more TSV's, which are the arrays of di/tri/tetranucleotide frequencies in each window. These files are large, especially as tetranucleotide frequencies will contain 4e4 columns. The kmers are sorted lexicographically from left -> right (AA(AA) to TT(TT)).

e.g. for dinucleotide frequencies:

```
ID	start	end	AA	AC	AG	AT	CA	CC	CG	CT	GA	GC	GG	GT	TA	TC	TG	TT
OV656674.1	0	1000	74	65	47	46	38	99	70	83	58	65	68	52	62	60	58	54
OV656674.1	1000	2000	36	54	58	33	49	89	120	49	75	101	97	65	22	63	62	26
OV656674.1	2000	3000	62	46	68	37	45	74	106	58	64	93	67	71	41	70	55	42
OV656674.1	3000	4000	55	56	55	59	61	59	81	58	70	73	72	66	39	71	73	51
OV656674.1	4000	5000	35	58	45	54	46	53	104	64	83	87	80	67	28	70	87	38
OV656674.1	5000	6000	81	71	42	40	29	106	68	87	56	60	71	51	68	53	58	58
OV656674.1	6000	7000	57	74	58	39	78	96	79	68	56	87	64	47	37	64	53	42
OV656674.1	7000	8000	81	45	79	52	63	74	65	62	69	82	59	52	44	62	59	51
OV656674.1	8000	9000	103	71	49	87	70	57	54	54	59	44	40	51	78	62	52	68
```

### Comments, updates & bugs

As of version 0.2.2, I've removed canonical kmers as an option; it was really computationally expensive and I couldn't think of a way to efficienty add it in. End users that wish this are pointed in the direction of <a href="https://github.com/tolkit/fw_group">fw_group</a>, which will at some point soon provide this functionality.

The masked (-m) flag only affects GC content, GC proportion, GC and AT skew, proportion of G's, C's, A's, T's, N's, CpG's. Kmers are coerced to uppercase automatically. Shannon index counts only uppercase nucleotides.

Please use, test, and let me know if there are any bugs or features you want implemented. Either raise an issue, or email me (see email in usage).
