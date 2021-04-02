pub mod utils {
    use std::str;

    pub const FOURFOLD_DEG: [&str; 32] = [
        "CTT", "CTA", "CTG", "CTC", "GTT", "GTC", "GTA", "GTG", "TCT", "TCC", "TCA", "TCG", "CCT",
        "CCC", "CCA", "CCG", "ACT", "ACC", "ACA", "ACG", "GCT", "GCC", "GCA", "GCG", "CGT", "CGC",
        "CGA", "CGG", "GGT", "GGC", "GGA", "GGG",
    ];

    pub fn gc_four(codons: Vec<&[u8]>) -> f64 {
        let codon_vec = codons
            .into_iter()
            .map(str::from_utf8)
            .collect::<Result<Vec<&str>, _>>()
            .unwrap_or(vec![]);

        let mut gc_counts = 0;

        for codon in codon_vec.iter() {
            let split_tuple = codon.split_at(2);
            if split_tuple.1 == "G" || split_tuple.1 == "C" {
                gc_counts += 1;
            }
        }
        let gc_four_float = gc_counts as f64 / codon_vec.len() as f64;
        gc_four_float
    }

    // from the frame trimming in main.rs, the sequence should start with full codons
    pub fn gc_four_fold_deg_sites(dna: &str) -> Vec<&[u8]> {
        //let four_fold_deg_sites = Vec::new();
        // convert to uppercase here, so bytes are the same
        let upper = dna.as_bytes();
        let codon_iterator = upper.chunks(3);
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
