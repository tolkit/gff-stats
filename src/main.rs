use bio::io::fasta;
use bio::io::gff;
use bio_types::strand::Strand;
use clap::{value_t, App, Arg};
use std::fs::{create_dir_all, File};
use std::io::{BufWriter, LineWriter, Write};
use std::str;
mod utils;
use rayon::prelude::*;
use std::sync::mpsc::channel;

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
        // .arg(
        //     Arg::with_name("sequences")
        //         .short("s")
        //         .long("sequences")
        //         .required(false)
        //         .default_value("true")
        //         .help("Save the extracted CDS fasta sequences?"),
        // )
        .get_matches();
    // parse command line options
    let input_gff = matches.value_of("gff").unwrap();
    let input_fasta = matches.value_of("fasta").unwrap();
    //let save_sequences = value_t!(matches.value_of("sequences"), bool).unwrap_or_else(|e| e.exit());

    let fasta_reader = fasta::Reader::from_file(input_fasta).expect("[-]\tPath invalid.");

    // create directory for output
    if let Err(e) = create_dir_all("./gff-stats/") {
        println!("[-]\tCreate directory error: {}", e.to_string());
    }
    let extension = "fasta";

    // create file
    let file_name = format!("./gff-stats/{}{}", "CDS.", extension);
    let gff_output = File::create(&file_name).expect("[-]\tUnable to create file");
    let gff_output = LineWriter::new(gff_output);
    let mut f = BufWriter::new(gff_output);

    // iterate over the fasta records
    // might have to do a manual loop.

    struct Output {
        fasta_id: String,
        start: usize,
        end: usize,
        gc_four: f64,
        seq: String,
    }

    // parallelise the processing of the fastas.
    let (sender, receiver) = channel();
    println!("[+]\tProcessing the following fastas in parallel:");
    fasta_reader
        .records()
        .par_bridge()
        .for_each_with(sender, |s, record| {
            let fasta_record = record.expect("[-]\tError during fasta record parsing.");
            println!(
                "[+]\t\t{} - {}bp.",
                fasta_record.id(),
                fasta_record.seq().len()
            );

            // now iterate over the gff file, matching on chromosome/scaff/contig ID
            let mut gff_reader =
                gff::Reader::from_file(input_gff, gff::GffType::GFF3).expect("[-]\tPath invalid.");

            // can we parallelise this too?
            for gff_record in gff_reader.records() {
                let rec = gff_record
                    .ok()
                    .expect("[-]\tError during gff record parsing.");

                if rec.seqname() != fasta_record.id() {
                    continue;
                } else if rec.seqname() == fasta_record.id() {
                    // only interested in coding sequences.
                    if rec.feature_type() != "CDS" {
                        continue;
                    } else if rec.feature_type() == "CDS" {
                        // trim the sequence, using start and end
                        // get start so the sequence always starts on the first base of a codon.
                        let start = *rec.start() as usize;
                        let end = *rec.end() as usize;
                        let slice = utils::utils::trim_sequence(
                            fasta_record.seq(),
                            start,
                            end,
                            rec.frame(),
                        );

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

                            s.send(Output {
                                fasta_id: fasta_record.id().to_string(),
                                start: start,
                                end: end,
                                gc_four: gc_four,
                                seq: revcomp_slice_str,
                            })
                            .expect("send!")
                        } else {
                            let four_deg_codons = utils::utils::gc_four_fold_deg_sites(&slice_str);
                            let gc_four = utils::utils::gc_four(four_deg_codons);

                            s.send(Output {
                                fasta_id: fasta_record.id().to_string(),
                                start: start,
                                end: end,
                                gc_four: gc_four,
                                seq: slice_str.to_owned(),
                            })
                            .expect("send!")
                        }
                    }
                }
            }
        });
    // collect and write.
    let res: Vec<Output> = receiver.iter().collect();

    println!("[+]\tCollected output and writing.");
    for i in res {
        writeln!(
            f,
            ">{}.{}-{}:{}\n{}",
            i.fasta_id, i.start, i.end, i.gc_four, i.seq
        )
        .unwrap_or_else(|_| println!("[-]\tError in writing to file."));
    }
    println!("[+]\tResults written.");
}
