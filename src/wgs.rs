pub mod wgs {

    pub struct GenomeStats {
        pub no_reads: usize,
        pub genome_length: usize,
        pub n50: usize, // collect N10 -> N50
    }
    // see https://github.com/rust-bio/rust-bio-tools/blob/e4ec0a88c6bf4bf804bdadbe58a30fd50571f2dd/src/sequences_stats.rs
    fn n_50(read_lengths: Vec<usize>, genome_length: usize) -> usize {
        let mut acc = 0;
        for val in read_lengths.iter() {
            acc += *val;
            if acc > genome_length / 2 {
                return *val;
            }
        }
        read_lengths[read_lengths.len() - 1]
    }

    pub fn collect_genome_stats(mut read_lengths: Vec<usize>) -> GenomeStats {
        let no_reads = read_lengths.len();
        let genome_length: usize = read_lengths.iter().sum();

        // calculate N50
        read_lengths.sort();
        let n50 = n_50(read_lengths, genome_length);
        // return the struct
        GenomeStats {
            no_reads: no_reads,
            genome_length: genome_length,
            n50: n50,
        }
    }
}
