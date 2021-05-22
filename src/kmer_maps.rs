pub mod kmer_maps {

    use std::collections::HashMap;
    use std::fmt::{Display, Error, Formatter};

    #[derive(Debug, Clone)]
    pub struct KmerMap {
        pub len: usize,
        pub map: HashMap<Vec<u8>, i32>,
    }

    pub fn generate_kmer_maps() -> Vec<KmerMap> {
        let kmer_maps: Vec<KmerMap> = vec![2, 3, 4]
            .iter()
            .map(|i| {
                let kmer_i = gen_all_kmers(*i);
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

    // display for Vec<i32>
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
}
