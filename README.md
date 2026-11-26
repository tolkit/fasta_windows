# fasta_windows

<p align="center">
    <img width="300" height="300" src="https://www.sanger.ac.uk/wp-content/uploads/DToL-logo_notext-300x300.png">
</p>

Written for Darwin Tree of Life chromosomal level genome assemblies. The executable takes a fasta formatted file and calculates some statistics of interest in windows:

- GC proportion
- GC and AT skew
- Proportion of G's, C's, A's, T's, N's, soft-masked bases, and CpG's
- Shannon entropy
- Context-Tree Weighting compressibility (optional, `-c`)
- Di/tri/tetranucleotide shannon diversity
- Mono/di/tri/tetranucleotide count arrays

An entropy-only mode (`-e`) is also available; see [below](#entropy-mode-centromere--repeat-detection-pipeline).

Output files can be visualised using <a href="https://github.com/tolkit/fw_plot">fw_plot</a> or grouped using <a href="https://github.com/tolkit/fw_group">fw_group</a>.

## Download

The easiest way to get `fasta_windows` is through <b><a href="https://anaconda.org/bioconda/fasta_windows">conda/bioconda</a></b>.

```bash
conda create -n fasta_windows -c bioconda fasta_windows
```

Or, with Rust installed, from <a href="https://crates.io/crates/fasta_windows">crates.io</a>:

```bash
cargo install fasta_windows
```

A macOS binary is also attached to each <a href="https://github.com/tolkit/fasta_windows/releases">GitHub release</a>. See [CHANGELOG.md](CHANGELOG.md) for what changed in each version.

## Usage

```
Quickly compute statistics over a fasta file in windows.

Usage: fasta_windows [OPTIONS] --fasta <fasta> --output <output>

Options:
  -f, --fasta <fasta>              The input fasta file.
  -w, --window_size <window_size>  Integer size of window for statistics to be computed over. [default: 1000]
  -d, --description                Add an extra column to _windows.tsv output with fasta header descriptions.
  -o, --output <output>            Output filename for the TSV's (without extension).
  -m, --masked                     Consider only uppercase nucleotides in the calculations.
  -c, --ctw                        Calculate the Context-Tree Weighting (how compressible a sequence is).
  -e, --entropy                    Entropy mode: output a single BED file of Shannon entropy and CTW per window. Skips all k-mer and nucleotide-composition computation for maximum speed.
  -h, --help                       Print help
  -V, --version                    Print version
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

The input may be plain or gzip-compressed FASTA. For bzip2, xz and zstd input, build with `cargo build --release --features compression` (this compiles the C compression libraries).

## Output

Output goes to `./fw_out/`. The main file, `<output>_freq_windows.tsv`, is a tsv with bed-like format in the first three columns. Rows are written in the same order as the records in the input FASTA. `-d` adds a `description` column after `ID`, and `-c` adds a `ctw` column after `Shannon_entropy`. For example, on `tests/data/small.fa`:

```
ID	start	end	GC_prop	GC_skew	AT_skew	Shannon_entropy	Prop_Gs	Prop_Cs	Prop_As	Prop_Ts	Prop_Ns	Prop_masked	CpG_prop	Dinucleotide_Shannon	Trinucleotide_Shannon	Tetranucleotide_Shannon
seq1	0	1000	0.493	0.014	-0.041	1.999	0.250	0.243	0.243	0.264	0.000	0.106	0.064	3.992	5.951	7.801
seq1	1000	2000	0.514	0.080	-0.017	2.307	0.236	0.201	0.203	0.210	0.150	0.263	0.052	3.989	5.957	7.797
seq1	2000	2500	0.478	-0.038	-0.080	1.996	0.230	0.248	0.240	0.282	0.000	0.334	0.058	3.980	5.907	7.594
seq2	0	700	0.501	-0.026	0.054	1.999	0.244	0.257	0.263	0.236	0.000	0.000	0.064	3.990	5.951	7.749
seq3	0	1000	0.300	1.000	0.429	1.485	0.300	0.000	0.500	0.200	0.000	0.000	0.000	2.331	2.541	2.550
```

Also output (non-optional at the moment) are four more TSVs, `<output>_{mononuc,dinuc,trinuc,tetranuc}_windows.tsv`, with the counts of each nucleotide (A, C, G, T, N) and each di/tri/tetranucleotide in each window. These files are large, especially the tetranucleotide file with 256 count columns. The kmers are sorted lexicographically from left -> right (AA(AA) to TT(TT)). Only k-mers made entirely of A, C, G and T (in either case) are counted; any other character, such as N or an IUPAC ambiguity code, breaks the k-mer.

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

## Entropy mode (centromere / repeat detection pipeline)

The `-e` flag enables a fast entropy-only mode that outputs a 5-column BED file suitable for repeat and centromere detection:

```bash
fasta_windows -e -f genome.fa -w 5000 -o species_name
# → fw_out/species_name_entropy.bed
# Columns: chrom  start  end  shannon_entropy  ctw_bits_per_base
```

Column 4 is Shannon entropy (bits, max 2.0 for equal ACGT). Column 5 is CTW (Context-Tree Weighting) bits per base — a compressibility measure that is independently sensitive to tandem repeats. Low values in either column indicate low-complexity sequence (satellite arrays, telomeres, centromeres).

This output feeds directly into **[fw_regions](https://github.com/tolkit/fw_regions)**, which detects the low-complexity regions and provides a satellite repeat database pipeline. See that repository for the full workflow.

A visualisation script is included:

```bash
Rscript plot_entropy.R -i fw_out/species_name_entropy.bed -o species_name
# Optional overlays: -r regions.bed -s satellite_windows.bed -m satellite_metadata.tsv
# Output: species_name_entropy_plot.pdf (one page per chromosome)
```

### Comments, updates & bugs

As of version 0.2.2, I've removed canonical kmers as an option; it was really computationally expensive and I couldn't think of a way to efficienty add it in. End users that wish this are pointed in the direction of <a href="https://github.com/tolkit/fw_group">fw_group</a>, which will at some point soon provide this functionality.

The masked (-m) flag makes GC proportion, GC and AT skew, the proportions of G's, C's, A's, T's and N's, and the mononucleotide counts consider only uppercase bases, and sets `Prop_masked` to 0. K-mers (and so CpG and the k-mer Shannon diversities) and `Shannon_entropy` ignore case either way. In entropy mode (`-e`), `-m` makes Shannon entropy count only uppercase bases.

Please use, test, and let me know if there are any bugs or features you want implemented. Either raise an issue, or email me (see `authors` in Cargo.toml).
