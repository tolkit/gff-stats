pub mod utils {

    use std::collections::HashMap;
    use std::str;

    pub const FOURFOLD_DEG: [&str; 32] = [
        "CTT", "CTA", "CTG", "CTC", "GTT", "GTC", "GTA", "GTG", "TCT", "TCC", "TCA", "TCG", "CCT",
        "CCC", "CCA", "CCG", "ACT", "ACC", "ACA", "ACG", "GCT", "GCC", "GCA", "GCG", "CGT", "CGC",
        "CGA", "CGG", "GGT", "GGC", "GGA", "GGG",
    ];

    // stats on the whole sequence
    pub struct Stats {
        pub gc_percent: f64,
        pub at_percent: f64,
        pub gc_skew: f64,
        pub at_skew: f64,
    }

    pub fn whole_seq_stats(dna: &str) -> Stats {
        // chars uses quite a bit more memory here
        // but sequences should never be super long...
        let map = dna.chars().fold(HashMap::new(), |mut map, c| {
            *map.entry(c).or_insert(0) += 1;
            map
        });
        let g_counts = map.get(&'G').unwrap_or(&0);
        let c_counts = map.get(&'C').unwrap_or(&0);
        let a_counts = map.get(&'A').unwrap_or(&0);
        let t_counts = map.get(&'T').unwrap_or(&0);

        Stats {
            gc_percent: (g_counts + c_counts) as f64 / dna.len() as f64,
            at_percent: (a_counts + t_counts) as f64 / dna.len() as f64,
            gc_skew: (g_counts - c_counts) as f64 / (g_counts + c_counts) as f64,
            at_skew: (a_counts - t_counts) as f64 / (a_counts + t_counts) as f64,
        }
    }

    // stats on the 4-fold-degenerate-sites
    pub struct FourFoldStats {
        pub gc_percent: f64,
        pub at_percent: f64,
        pub gc_skew: f64,
        pub at_skew: f64,
    }
    pub fn four_fold_site_stats(codons: Vec<&[u8]>) -> FourFoldStats {
        let codon_vec = codons
            .into_iter()
            .map(str::from_utf8)
            .collect::<Result<Vec<&str>, _>>()
            .unwrap_or(vec![]);

        let mut g_counts = 0;
        let mut c_counts = 0;
        let mut a_counts = 0;
        let mut t_counts = 0;

        for codon in codon_vec.iter() {
            let split_tuple = codon.split_at(2);

            match split_tuple.1 {
                "G" => g_counts += 1,
                "C" => c_counts += 1,
                "A" => a_counts += 1,
                "T" => t_counts += 1,
                _ => (),
            }
        }
        FourFoldStats {
            gc_percent: (g_counts + c_counts) as f64 / codon_vec.len() as f64,
            at_percent: (a_counts + t_counts) as f64 / codon_vec.len() as f64,
            gc_skew: (g_counts - c_counts) as f64 / (g_counts + c_counts) as f64,
            at_skew: (a_counts - t_counts) as f64 / (a_counts + t_counts) as f64,
        }
    }

    // from the frame trimming in main.rs, the sequence should start with full codons
    pub fn gc_four_fold_deg_sites(dna: &str) -> Vec<&[u8]> {
        let bytes = dna.as_bytes();
        let codon_iterator = bytes.chunks(3);
        let mut collector = Vec::new();

        for codon in codon_iterator {
            let str_codon = str::from_utf8(codon).unwrap_or("");
            if FOURFOLD_DEG.iter().any(|&i| i == str_codon.to_uppercase()) {
                collector.push(codon);
            }
        }
        collector
    }

    pub fn trim_sequence<'a>(seq: &'a [u8], start: usize, end: usize, frame: &str) -> &'a [u8] {
        let mut trimmed: &'a [u8] = seq.get(start..end).unwrap_or(b"");

        if frame == "1" {
            trimmed = trimmed.get(2..).unwrap_or(b"");
        } else if frame == "2" {
            trimmed = trimmed.get(1..).unwrap_or(b"");
        } else {
            trimmed = trimmed;
        }

        if trimmed.len() % 3 == 0 {
            return trimmed;
        } else if trimmed.len() % 3 == 1 {
            return trimmed.get(..trimmed.len() - 1).unwrap_or(b"");
        } else {
            return trimmed.get(..trimmed.len() - 2).unwrap_or(b"");
        }
    }

    pub fn reverse_complement(dna: &str) -> String {
        let dna_chars = dna.chars();
        let mut revcomp = Vec::new();

        for base in dna_chars {
            revcomp.push(switch_base(base))
        }
        revcomp.as_mut_slice().reverse();
        revcomp.into_iter().collect()
    }

    // switch lowercase letters too here.
    fn switch_base(c: char) -> char {
        match c {
            'A' => 'T',
            'a' => 't',
            'C' => 'G',
            'c' => 'g',
            'T' => 'A',
            't' => 'a',
            'G' => 'C',
            'g' => 'c',
            'N' => 'N',
            'n' => 'n',
            _ => 'N',
        }
    }
}
