use bio::io::fasta;
use bio::io::gff;
use bio_types::strand::Strand;
use clap::{App, Arg};
use std::str;

mod utils;

fn main() {
    // command line options
    let matches = App::new("Fasta windows")
        .version(clap::crate_version!())
        .author("Max Brown <mb39@sanger.ac.uk>")
        .about("Extract GFF3 regions and compute statistics on them.")
        .arg(
            Arg::with_name("gff")
                .short("g")
                .long("gff")
                .takes_value(true)
                .required(true)
                .help("The input gff file."),
        )
        .arg(
            Arg::with_name("fasta")
                .short("f")
                .long("fasta")
                .takes_value(true)
                .required(true)
                .help("The reference fasta file."),
        )
        .get_matches();
    // parse command line options
    let input_gff = matches.value_of("gff").unwrap();
    let input_fasta = matches.value_of("fasta").unwrap();

    // make reverse revcomp in fasta?

    let mut gff_reader =
        gff::Reader::from_file(input_gff, gff::GffType::GFF3).expect("[-]\tPath invalid.");
    let fasta_reader = fasta::Reader::from_file(input_fasta).expect("[-]\tPath invalid.");

    println!("[+]\tExtracting fasta files.");
    // iterate over the fasta records
    for result in fasta_reader.records() {
        let fasta_record = result.expect("[-]\tError during fasta record parsing.");
        // now iterate over the gff file, matching on chromosome/scaff/contig ID
        for gff_record in gff_reader.records() {
            let rec = gff_record
                .ok()
                .expect("[-]\tError during gff record parsing.");

            if rec.seqname() == fasta_record.id() {
                // only interested in coding sequences.
                if rec.feature_type() == "CDS" {
                    // trim the sequence, using start and end
                    // get start so the sequence always starts on the first base of a codon.
                    let start = *rec.start() as usize;
                    let end = *rec.end() as usize;
                    let slice =
                        utils::utils::trim_sequence(fasta_record.seq(), start, end, rec.frame());
                    //let slice = fasta_record.seq().get(start..end).unwrap();

                    let slice_str = str::from_utf8(slice).unwrap();
                    // if reverse, take the revcomp
                    let strandedness = match rec.strand() {
                        Some(strand) => strand,
                        None => Strand::Unknown,
                    };

                    if strandedness == Strand::Reverse {
                        let revcomp_slice_str = utils::utils::reverse_complement(slice_str);
                        let four_deg_codons =
                            utils::utils::gc_four_fold_deg_sites(&revcomp_slice_str);
                        let gc_four = utils::utils::gc_four(four_deg_codons);

                        // we can save the intermediate fastas if we want?
                        println!(
                            ">{}.{}-{}:{}\n{:?}",
                            fasta_record.id(),
                            start,
                            end,
                            gc_four,
                            revcomp_slice_str
                        );
                    } else {
                        let four_deg_codons = utils::utils::gc_four_fold_deg_sites(&slice_str);
                        let gc_four = utils::utils::gc_four(four_deg_codons);
                        // we can save the intermediate fastas if we want?
                        println!(
                            ">{}.{}-{}:{}\n{:?}",
                            fasta_record.id(),
                            start,
                            end,
                            gc_four,
                            slice_str
                        );
                    }
                }
            }
        }
    }
}
