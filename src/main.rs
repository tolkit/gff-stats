use bio::io::fasta;
use bio::io::gff;
use bio_types::strand::Strand;
use clap::{value_t, App, Arg};
use std::fs::{create_dir_all, File};
use std::io::LineWriter;
use std::str;
mod utils;

fn main() {
    // command line options
    let matches = App::new("GFF stats")
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
        .arg(
            Arg::with_name("sequences")
                .short("s")
                .long("sequences")
                .required(false)
                .default_value("true")
                .help("Save the extracted CDS fasta sequences?"),
        )
        .get_matches();
    // parse command line options
    let input_gff = matches.value_of("gff").unwrap();
    let input_fasta = matches.value_of("fasta").unwrap();
    let save_sequences = value_t!(matches.value_of("sequences"), bool).unwrap_or_else(|e| e.exit());

    let fasta_reader = fasta::Reader::from_file(input_fasta).expect("[-]\tPath invalid.");

    // create directory for output
    if let Err(e) = create_dir_all("./gff-stats/") {
        println!("[-]\tCreate directory error: {}", e.to_string());
    }
    let extension = "fasta";

    // create file
    let file_name = format!("./gff-stats/{}{}", "CDS.", extension);
    let gff_output = File::create(&file_name).unwrap();
    let mut gff_output = LineWriter::new(gff_output);

    println!("[+]\tExtracting fasta files.");
    // iterate over the fasta records
    // might have to do a manual loop.
    for result in fasta_reader.records() {
        let fasta_record = result.expect("[-]\tError during fasta record parsing.");
        println!(
            "[+]\tProcessing {} - {}bp.",
            fasta_record.id(),
            fasta_record.seq().len()
        );
        // now iterate over the gff file, matching on chromosome/scaff/contig ID
        let mut gff_reader =
            gff::Reader::from_file(input_gff, gff::GffType::GFF3).expect("[-]\tPath invalid.");
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
                        if save_sequences {
                            utils::utils::create_dir_file(
                                &mut gff_output,
                                fasta_record.id(),
                                start,
                                end,
                                gc_four,
                                &revcomp_slice_str,
                            )
                        }
                    } else {
                        let four_deg_codons = utils::utils::gc_four_fold_deg_sites(&slice_str);
                        let gc_four = utils::utils::gc_four(four_deg_codons);
                        // we can save the intermediate fastas if we want?
                        if save_sequences {
                            utils::utils::create_dir_file(
                                &mut gff_output,
                                fasta_record.id(),
                                start,
                                end,
                                gc_four,
                                slice_str,
                            )
                        }
                    }
                }
            } else {
                continue;
            }
        }
    }
}
