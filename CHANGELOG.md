# Changelog

All notable changes to fasta_windows are documented here. Versions follow
[semantic versioning](https://semver.org/).

## [0.3.0] - 2026-11-25 "Bothriochloa"

A speed release. Output is unchanged, except for the row order and the bug fixes
listed below: on inputs without IUPAC ambiguity codes, every output file is
byte-identical to 0.2.x once rows are sorted by ID.

Version 0.2.5 and 0.2.6 were never released. Their changes (CTW, entropy mode,
the masked proportion) ship for the first time in this release.

### Performance

Measured on an Apple M1 (8 threads) with synthetic genomes:

| Input        | Mode    | 0.2.6              | 0.3.0              | Speed-up |
| ------------ | ------- | ------------------ | ------------------ | -------- |
| 100 Mb       | default | 6.9 s              | 0.31 s             | 22×      |
| 100 Mb       | `-c`    | 25.8 s             | 0.93 s             | 28×      |
| 100 Mb       | `-e`    | 9.4 s              | 0.92 s             | 10×      |
| 1 Gb         | default | 43 s, 1.65 GB peak | 1.6 s, 0.43 GB peak | 27×     |

- k-mers are counted in a single pass with a rolling 2-bit code into fixed
  arrays, instead of allocating and hashing a `Vec<u8>` per k-mer (~40×).
- CTW now only updates node counts per base and evaluates the tree once per
  window, instead of computing every mixture after every base (~12×). Results
  match the previous implementation to within 1e-12 bits per base.
- Nucleotide statistics come from one 256-bin count table instead of a
  per-window `HashMap` (~15×).
- Windows are processed in parallel within and across records, so a genome
  dominated by a few large chromosomes now uses every core.
- Input is parsed once (previously twice) on a reader thread, output is
  formatted in parallel and written on a writer thread, and memory stays
  bounded by batch size instead of growing with genome size.
- Release builds use fat LTO and a single codegen unit.

### Added

- Entropy mode, `-e/--entropy`: writes a BED file of Shannon entropy and CTW
  per window, for repeat and centromere detection. It is now built by default
  (the `entropy` Cargo feature is on by default). `plot_entropy.R` plots the
  output.
- `-c/--ctw`: a Context-Tree Weighting (compressibility) column in
  `_freq_windows.tsv`.
- `Prop_masked` column (proportion of soft-masked bases), and the W and S
  ambiguity codes count towards GC proportion (#12, Rich Challis).
- gzip-compressed FASTA input. bzip2, xz and zstd input need the new
  `compression` Cargo feature, which builds C libraries.

### Changed

- **Rows are written in FASTA order.** Previously all rows were sorted by
  record ID (`chr1, chr10, chr2, …`).
- FASTA parsing uses needletail instead of bio. The progress bar now shows
  bytes read rather than records.
- The minimum supported Rust version is 1.82.
- Library API: `kmer_maps` is removed; `kmeru8::count_kmers`,
  `kmeru8::shannon` and the new `pipeline` module replace it, and
  `SeqStats::nuc_counts` is now `[u64; 5]`.

### Fixed

- IUPAC ambiguity codes no longer corrupt the k-mer output. A k-mer containing
  any base other than A, C, G or T (in either case) is now skipped; previously
  only N was skipped, and k-mers like `AR` were added as extra columns. That
  misaligned those rows against the header, and `CpG_prop` read the wrong
  column.
- k-mer Shannon diversity was summed in `HashMap` order, so its last bits
  varied between runs.
- `cargo test` passes again; there are new golden-output tests that pin the
  output to 0.2.x.

### Removed

- The `bio` and `itertools` dependencies.

## [0.2.4] - 2022-08-23 "Laplacea"

- Mono-nucleotide counts, N counts, CpG proportion and AT skew in the output
  (Matthieu Muffato).
- Better error handling; updated indicatif.

## [0.2.3] - 2022-04-25 "Monolophus"

- macOS binary attached to the release.
- Updated clap, removed redundant modules, added unit tests.
- Fixed a windowing bug (Matthieu Muffato).

## [0.2.2] - 2021-11-02 "Conda release"

- First release on bioconda.
- Canonical k-mers removed (too computationally expensive); see
  [fw_group](https://github.com/tolkit/fw_group).

[0.3.0]: https://github.com/tolkit/fasta_windows/compare/v0.2.4...v0.3.0
[0.2.4]: https://github.com/tolkit/fasta_windows/releases/tag/v0.2.4
[0.2.3]: https://github.com/tolkit/fasta_windows/releases/tag/v0.2.3
[0.2.2]: https://github.com/tolkit/fasta_windows/releases/tag/v0.2.2
