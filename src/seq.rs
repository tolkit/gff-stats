pub mod seq {

    use bio::io::fasta;
    use bio::io::gff;
    use bio_types::strand::Strand;
    use rayon::prelude::*;
    use std::str;

    use crate::utils;

    // as far as I know, the printing of fastas is not repeatable, or fully ordered
    // as it's run in parallel

    pub fn generate_seqs(matches: &clap::ArgMatches) {
        // parse command line options
        let input_gff = matches.value_of("gff").unwrap();
        let input_fasta = matches.value_of("fasta").unwrap();

        // does the user want to print out spliced CDS?
        let spliced = matches.is_present("spliced");

        if spliced {
            run_spliced(input_gff, input_fasta)
        } else {
            run_cds(input_gff, input_fasta)
        }
    }

    // much faster.
    fn run_cds(input_gff: &str, input_fasta: &str) {
        let fasta_reader = fasta::Reader::from_file(input_fasta).expect("[-]\tPath invalid.");

        eprintln!("[+]\tProcessing the following fasta sequences in parallel:");
        fasta_reader.records().par_bridge().for_each(|record| {
            let fasta_record = record.expect("[-]\tError during fasta record parsing.");
            eprintln!(
                "[+]\t\t{} - {}bp",
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

                        let attr = rec.attributes();

                        let slice = utils::utils::trim_sequence(
                            fasta_record.seq(),
                            start,
                            end,
                            rec.frame(),
                        );

                        // TODO: check this is necessary? Will require changing
                        // utils functions.
                        let slice_str = str::from_utf8(slice).unwrap();
                        // if reverse, take the revcomp
                        let strandedness = match rec.strand() {
                            Some(strand) => strand,
                            None => Strand::Unknown,
                        };

                        if strandedness == Strand::Reverse {
                            let revcomp_slice_str = &utils::utils::reverse_complement(slice_str);

                            println!(
                                ">{}: {}(-): {}-{}\n{}",
                                rec.seqname(),
                                attr.get("ID").unwrap(),
                                start,
                                end,
                                revcomp_slice_str
                            )
                        } else {
                            println!(
                                ">{}: {}(+): {}-{}\n{}",
                                rec.seqname(),
                                attr.get("ID").unwrap(),
                                start,
                                end,
                                slice_str
                            )
                        }
                    }
                }
            }
        });
        eprintln!("[+]\tFinished.");
    }

    fn run_spliced(input_gff: &str, input_fasta: &str) {
        let fasta_reader = fasta::Reader::from_file(input_fasta).expect("[-]\tPath invalid.");
        #[derive(Debug)]
        struct SplicedCDS {
            seq: String,
            start: usize,
            end: usize,
        }

        eprintln!("[+]\tProcessing the following fasta sequences in parallel:");
        fasta_reader.records().par_bridge().for_each(|record| {
            let fasta_record = record.expect("[-]\tError during fasta record parsing.");
            eprintln!(
                "[+]\t\t{} - {}bp",
                fasta_record.id(),
                fasta_record.seq().len()
            );

            // now iterate over the gff file, matching on chromosome/scaff/contig ID
            let mut gff_reader_1 =
                gff::Reader::from_file(input_gff, gff::GffType::GFF3).expect("[-]\tPath invalid.");
            let mut gff_reader_2 =
                gff::Reader::from_file(input_gff, gff::GffType::GFF3).expect("[-]\tPath invalid.");

            // we need to look ahead to the next iteration.
            // push
            let mut spliced: Vec<SplicedCDS> = Vec::new();

            for (gff_record_1, gff_record_2) in
                gff_reader_1.records().zip(gff_reader_2.records().skip(1))
            {
                let rec_1 = gff_record_1
                    .ok()
                    .expect("[-]\tError during gff record parsing.");
                let rec_2 = gff_record_2
                    .ok()
                    .expect("[-]\tError during gff record parsing.");

                if rec_1.seqname() != fasta_record.id() {
                    continue;
                } else if rec_1.seqname() == fasta_record.id() {
                    // only interested in coding sequences.
                    if rec_1.feature_type() != "CDS" {
                        continue;
                    } else if rec_1.feature_type() == "CDS" {
                        // attributes
                        let attr_1 = rec_1.attributes();
                        let attr_2 = rec_2.attributes();

                        // trim the sequence, using start and end
                        // get start so the sequence always starts on the first base of a codon.
                        let start = *rec_1.start() as usize;
                        let end = *rec_1.end() as usize;

                        let slice_1 = utils::utils::trim_sequence(
                            fasta_record.seq(),
                            start,
                            end,
                            rec_1.frame(),
                        );
                        let slice_2 = utils::utils::trim_sequence(
                            fasta_record.seq(),
                            start,
                            end,
                            rec_1.frame(),
                        );

                        // TODO: check this is necessary? Will require changing
                        // utils functions.
                        let slice_str_1 = str::from_utf8(slice_1).unwrap();
                        let slice_str_2 = str::from_utf8(slice_2).unwrap();
                        // if reverse, take the revcomp
                        let strandedness = match rec_1.strand() {
                            Some(strand) => strand,
                            None => Strand::Unknown,
                        };

                        if strandedness == Strand::Reverse {
                            let revcomp_slice_str_1 =
                                &utils::utils::reverse_complement(slice_str_1);
                            let revcomp_slice_str_2 = utils::utils::reverse_complement(slice_str_2);

                            if attr_1.get("ID").unwrap() == attr_2.get("ID").unwrap() {
                                spliced.push(SplicedCDS {
                                    seq: revcomp_slice_str_1.clone(),
                                    start,
                                    end,
                                })
                            } else {
                                // add the final
                                spliced.push(SplicedCDS {
                                    seq: revcomp_slice_str_2.clone(),
                                    start,
                                    end,
                                });

                                use std::cmp::Reverse;
                                spliced.sort_unstable_by_key(|d| Reverse(d.start));

                                let seq: Vec<&str> =
                                    spliced.iter().map(|x| x.seq.as_str()).collect::<Vec<_>>();
                                let positions: Vec<String> = spliced
                                    .iter()
                                    .map(|x| format!("{}-{}", x.start, x.end))
                                    .collect::<Vec<_>>();
                                let res = seq.join("");
                                let res2 = positions.join(",");

                                println!(
                                    ">{}: {}(-): {}\n{}",
                                    rec_1.seqname(),
                                    attr_1.get("ID").unwrap(),
                                    res2,
                                    res
                                );
                                spliced.clear();
                            }
                        } else {
                            // if strandedness is forward
                            if attr_1.get("ID").unwrap() == attr_2.get("ID").unwrap() {
                                spliced.push(SplicedCDS {
                                    seq: slice_str_1.to_string(),
                                    start,
                                    end,
                                })
                            } else {
                                // add the final
                                spliced.push(SplicedCDS {
                                    seq: slice_str_2.to_string(),
                                    start,
                                    end,
                                });

                                spliced.sort_unstable_by_key(|d| d.start);

                                let seq: Vec<&str> =
                                    spliced.iter().map(|x| x.seq.as_str()).collect::<Vec<_>>();
                                let positions: Vec<String> = spliced
                                    .iter()
                                    .map(|x| format!("{}-{}", x.start, x.end))
                                    .collect::<Vec<_>>();
                                let res = seq.join("");
                                let res2 = positions.join(",");

                                println!(
                                    ">{}: {}(+): {}\n{}",
                                    rec_1.seqname(),
                                    attr_1.get("ID").unwrap(),
                                    res2,
                                    res
                                );
                                spliced.clear();
                            }
                        }
                    }
                }
            }
        });
        eprintln!("[+]\tFinished.");
    }
}
