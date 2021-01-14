pub mod utils {

    // functions to aid kmer distance to reference in windows
    // and to wite kmer count spectra in windows
    
    use std::collections::HashMap;

    // relatively inflexible function which takes a reference hashmap and adds the values of the second,
    // or adds the key and value if no key is present in the reference, *in place*.

    pub fn merge_hashmap_ip(a: &mut HashMap<Vec<u8>, i32>, b: HashMap<Vec<u8>, i32>) -> &HashMap<Vec<u8>, i32> {
        // can I implement this without cloning here..? Probably a way of parameterising this differently
        // b will always be the current record though...
        let b_clone = b.clone();
        for (k, v) in b.into_iter() {
            if a.get(&k).is_none() {
                a.insert(k, v);
            } else {
                let count_to_add = b_clone.get(&k).unwrap();
                // update the value of the key
                *a.get_mut(&k).unwrap() += count_to_add;
            }
        }
        a
    }

    // a.k.a the kmer frequency count.
    fn freq(entries: Vec<i32>) -> HashMap<i32, i32> {
        entries.iter().fold(HashMap::new(), |mut freqs, value| {
        *freqs.entry(*value).or_insert(0) += 1;
        freqs
    })
    }

    // then we want to compute the euclidean distance
    // where vec1 is the window derived 
    fn euclidean_distance(v1: &Vec<i32>, v2: &Vec<i32>) -> f64 {
        let x: Vec<i32> = v1.iter()
            .zip(v2.iter())
            .map(|(x,y)| (x - y).pow(2)).collect();
        
        let y: f64 = x.iter().sum::<i32>().into();
        let res = y.sqrt();
        res
    }

    // input is the hashmaps of the entire genome of kmers & window kmers
    pub fn create_kmer_distance(window: HashMap<Vec<u8>, i32>, reference: &mut HashMap<Vec<u8>, i32>) -> f64 {
        // get values as vectors
        let window_values = window.values().cloned().collect::<Vec<i32>>();
        let reference_values = reference.values().cloned().collect::<Vec<i32>>();
        // tabulate kmer vector
        let window_kmer_table = freq(window_values);
        let reference_kmer_table = freq(reference_values);
        // finally extract the values.
        let window_kmer_table_values = window_kmer_table.values().cloned().collect::<Vec<i32>>();
        let reference_kmer_table_values = reference_kmer_table.values().cloned().collect::<Vec<i32>>();

        euclidean_distance(&window_kmer_table_values, &reference_kmer_table_values)
    }
}
