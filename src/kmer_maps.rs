pub mod kmer_maps {

    use std::collections::HashMap;
    use std::fmt::{Display, Error, Formatter};

    use crate::kmeru8::kmeru8::reverse_complement;

    #[derive(Debug, Clone)]
    pub struct KmerMap {
        pub len: usize,
        pub map: HashMap<Vec<u8>, i32>,
    }

    pub fn generate_kmer_maps(canonical: bool) -> Vec<KmerMap> {
        let kmer_maps: Vec<KmerMap> = vec![2, 3, 4]
            .iter()
            .map(|i| {
                let mut kmer_i = gen_all_kmers(*i);

                // if canonical = true, call filter_canonical
                match canonical {
                    true => {
                        kmer_i = filter_canonical(kmer_i);
                    }
                    false => (),
                }

                let mut kmers_u8: Vec<Vec<u8>> = Vec::new();
                for i in kmer_i {
                    kmers_u8.push(i.as_bytes().to_vec());
                }
                let mut map = HashMap::new();
                // insert all the possible kmers as keys
                for k in kmers_u8 {
                    map.insert(k, 0i32);
                }

                KmerMap { len: *i, map }
            })
            .collect();
        kmer_maps
    }

    // filter Vec<String> of kmers for canonical only.
    // so the hashmap will only contain canonical kmer keys

    fn filter_canonical(kmers: Vec<String>) -> Vec<String> {
        let revcomp_kmers: Vec<String> = kmers
            .iter()
            .map(|e| {
                std::str::from_utf8(&reverse_complement(e.as_bytes()))
                    .unwrap()
                    .to_owned()
            })
            .collect();

        let mut canonical_kmers = Vec::new();

        for (e1, e2) in kmers.into_iter().zip(revcomp_kmers.into_iter()) {
            if e1 < e2 {
                canonical_kmers.push(e1);
            } else {
                canonical_kmers.push(e2);
            }
        }
        canonical_kmers.sort_unstable();
        canonical_kmers.dedup();
        canonical_kmers
    }

    // move out of kmeru8.rs
    fn gen_all_kmers(k: usize) -> Vec<String> {
        let mut output: Vec<String> = Vec::new();
        let nucleotides = ["A", "C", "G", "T"].iter().map(|x| x.to_string()).collect();
        if k == 0 {
            return output;
        } else if k == 1 {
            return nucleotides;
        } else {
            for suffix in gen_all_kmers(k - 1) {
                for nt in &nucleotides {
                    let new = suffix.clone() + &nt;
                    output.push(new);
                }
            }
        }
        output
    }

    // hacky display for Vec<i32>
    #[derive(Clone)]
    pub struct WriteArray(pub Vec<i32>);

    impl Display for WriteArray {
        fn fmt(&self, f: &mut Formatter) -> Result<(), Error> {
            let mut tab_separated = String::new();

            for num in &self.0[0..self.0.len() - 1] {
                tab_separated.push_str(&num.to_string());
                tab_separated.push_str("\t");
            }

            tab_separated.push_str(&self.0[self.0.len() - 1].to_string());
            write!(f, "{}", tab_separated)
        }
    }
    // hacky display for Vec<Vec<u8>>
    // which is what the kmers are stored as for most of the time.
    #[derive(Clone)]
    pub struct WriteKmerValues<'a>(pub Vec<&'a Vec<u8>>);

    impl<'a> Display for WriteKmerValues<'a> {
        fn fmt(&self, f: &mut Formatter) -> Result<(), Error> {
            let mut tab_separated = String::new();

            for kmer in &self.0[0..self.0.len() - 1] {
                let kmer_str = std::str::from_utf8(kmer).unwrap();
                tab_separated.push_str(kmer_str);
                tab_separated.push_str("\t");
            }

            tab_separated.push_str(std::str::from_utf8(&self.0[self.0.len() - 1]).unwrap());
            write!(f, "{}", tab_separated)
        }
    }
}
