use anyhow::Result;
use bio::io::fasta;
use bio::io::gff;
use bio_types::strand::Strand;
use rayon::prelude::*;
use std::cmp::Reverse;
use std::fs::{create_dir_all, File};
use std::io::{BufWriter, LineWriter, Write};
use std::str;
use std::sync::mpsc::channel;

use crate::utils;

/// Wrapper function to calculate statistics on sequences from a
/// GFF3 file.
pub fn calculate_stats(matches: &clap::ArgMatches) -> Result<()> {
    // parse command line options
    let input_gff = matches.value_of("gff").unwrap();
    let input_fasta = matches.value_of("fasta").unwrap();
    let degeneracy = matches.value_of_t("degeneracy").unwrap();
    let spliced = matches.is_present("spliced");
    let output = matches.value_of("output").unwrap();

    // create directory for output
    if let Err(e) = create_dir_all("./gff-stats/") {
        eprintln!("[-]\tCreate directory error: {}", e.to_string());
    }

    let file_name = format!("./gff-stats/{}_stats.tsv", output);
    let stats = File::create(&file_name)?;
    let stats = LineWriter::new(stats);
    let mut f = BufWriter::new(stats);

    if spliced {
        calculate_stats_spliced(input_gff, input_fasta, degeneracy, &mut f)?
    } else {
        calculate_stats_unspliced(input_gff, input_fasta, degeneracy, &mut f)?
    }

    Ok(())
}

/// Data structure sent through the parallel iterator.
pub struct Output {
    /// ID of the fasta
    pub fasta_id: String,
    /// ID of the GFF3 attribute
    pub attr_id: String,
    /// Index start of the sequence
    pub start: usize,
    /// Index end of the sequence
    pub end: usize,
    /// The statistics generated on the sequence
    pub stats: utils::Stats,
    /// The four-fold-degenerate statistics calculated on the sequences
    pub gc_four_stats: utils::FourFoldStats,
    /// The GC content of the third codon position
    pub gc_3_stats: utils::GC3Stats,
}

/// Wrapper for a vector of type [`Output`].
pub struct Writer {
    output: Vec<Output>,
}

/// The headers of the output table.
const HEADERS: [&str; 16] = [
    "fasta_id",
    "attr_id",
    "start",
    "end",
    "four_fold_deg_gc_per",
    "four_fold_deg_at_per",
    "four_fold_deg_gc_skew",
    "four_fold_deg_at_skew",
    "gc_per",
    "at_per",
    "gc_skew",
    "at_skew",
    "gc3_per",
    "at3_per",
    "gc3_sker",
    "at3_skew",
];

impl Writer {
    /// Write a vector of type [`Output`].
    pub fn write<T: Write>(&self, file: &mut BufWriter<LineWriter<T>>) -> Result<()> {
        let headers = HEADERS.join("\t");
        // write the headers
        writeln!(file, "{}", headers)?;

        for i in &self.output {
            writeln!(
                file,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                i.fasta_id,
                i.attr_id,
                i.start,
                i.end,
                i.gc_four_stats.gc_percent,
                i.gc_four_stats.at_percent,
                i.gc_four_stats.gc_skew,
                i.gc_four_stats.at_skew,
                i.stats.gc_percent,
                i.stats.at_percent,
                i.stats.gc_skew,
                i.stats.at_skew,
                i.gc_3_stats.gc_percent,
                i.gc_3_stats.at_percent,
                i.gc_3_stats.gc_skew,
                i.gc_3_stats.at_skew
            )?;
        }
        eprintln!("[+]\tResults written.");

        Ok(())
    }
}

/// Non-spliced version of the function which prints out statistics to
/// `STDOUT`.
///
/// Note the order of output statistics is not repeatable, as they are
/// processed in parallel.
fn calculate_stats_unspliced<T: Write>(
    input_gff: &str,
    input_fasta: &str,
    degeneracy: String,
    f: &mut BufWriter<LineWriter<T>>,
) -> Result<()> {
    let fasta_reader = fasta::Reader::from_file(input_fasta)?;

    // parallelise the processing of the fastas.
    let (sender, receiver) = channel();
    eprintln!("[+]\tProcessing the following fasta sequences in parallel:");

    fasta_reader
        .records()
        .par_bridge()
        .for_each_with(sender, |s, record| {
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
                        // attributes
                        let attr = rec.attributes();

                        // trim the sequence, using start and end
                        // get start so the sequence always starts on the first base of a codon.
                        let start = *rec.start() as usize;
                        let end = *rec.end() as usize;

                        let slice = utils::trim_sequence(
                            fasta_record.seq(),
                            start,
                            end,
                            rec.frame(),
                            false,
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
                            let revcomp_slice_str = &utils::reverse_complement(slice_str);

                            // I presume doing stats on the spliced seq would go here.
                            let gc_3_stats = utils::gc_3(
                                revcomp_slice_str,
                                Some(attr.get("ID").unwrap().clone()),
                            );
                            let seq_stats = utils::whole_seq_stats(revcomp_slice_str);
                            let four_deg_codons =
                                utils::gc_four_fold_deg_sites(revcomp_slice_str, &degeneracy);
                            let gc_four = utils::four_fold_site_stats(four_deg_codons);

                            s.send(Output {
                                fasta_id: fasta_record.id().to_string(),
                                attr_id: attr.get("ID").unwrap().clone(),
                                start,
                                end,
                                stats: seq_stats,
                                gc_four_stats: gc_four,
                                gc_3_stats,
                            })
                            .expect("Did not send.");
                        } else {
                            let gc_3_stats =
                                utils::gc_3(slice_str, Some(attr.get("ID").unwrap().clone()));
                            let seq_stats = utils::whole_seq_stats(slice_str);
                            let four_deg_codons =
                                utils::gc_four_fold_deg_sites(slice_str, &degeneracy);
                            let gc_four = utils::four_fold_site_stats(four_deg_codons);

                            s.send(Output {
                                fasta_id: fasta_record.id().to_string(),
                                attr_id: attr.get("ID").unwrap().clone(),
                                start,
                                end,
                                stats: seq_stats,
                                gc_four_stats: gc_four,
                                gc_3_stats,
                            })
                            .expect("Did not send.");
                        }
                    }
                }
            }
        });
    // collect and write.
    let res: Vec<Output> = receiver.iter().collect();
    let writer = Writer { output: res };

    eprintln!("[+]\tCollected output and writing.");

    writer.write(f)?;

    Ok(())
}

/// Store the spliced CDS sequences
#[derive(Debug)]
pub struct SplicedCDS {
    /// The sequence of the CDS
    pub seq: String,
    /// Index start of the sequence
    pub start: usize,
    /// Index end of the sequence
    pub end: usize,
}

/// Spliced version of the function which prints out statistics to
/// `STDOUT`.
///
/// Note the order of output statistics is not repeatable, as they are
/// processed in parallel.
fn calculate_stats_spliced<T: Write>(
    input_gff: &str,
    input_fasta: &str,
    degeneracy: String,
    f: &mut BufWriter<LineWriter<T>>,
) -> Result<()> {
    let fasta_reader = fasta::Reader::from_file(input_fasta)?;

    // parallelise the processing of the fastas.
    let (sender, receiver) = channel();
    eprintln!("[+]\tProcessing the following fasta sequences in parallel:");

    fasta_reader
        .records()
        .par_bridge()
        .for_each_with(sender, |s, record| {
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

                        let slice_1 = utils::trim_sequence(
                            fasta_record.seq(),
                            start,
                            end,
                            rec_1.frame(),
                            true,
                        );
                        let slice_2 = utils::trim_sequence(
                            fasta_record.seq(),
                            start,
                            end,
                            rec_1.frame(),
                            true,
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
                            let revcomp_slice_str_1 = &utils::reverse_complement(slice_str_1);
                            let revcomp_slice_str_2 = utils::reverse_complement(slice_str_2);

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

                                spliced.sort_unstable_by_key(|d| Reverse(d.start));

                                let seq: Vec<&str> =
                                    spliced.iter().map(|x| x.seq.as_str()).collect::<Vec<_>>();
                                let min_start: usize = *spliced
                                    .iter()
                                    .map(|x| x.start)
                                    .collect::<Vec<_>>()
                                    .first()
                                    .unwrap();
                                let max_end: usize = *spliced
                                    .iter()
                                    .map(|x| x.end)
                                    .collect::<Vec<_>>()
                                    .last()
                                    .unwrap();
                                let res = seq.join("");

                                // I presume doing stats on the spliced seq would go here.
                                let gc_3_stats =
                                    utils::gc_3(&res, Some(attr_1.get("ID").unwrap().clone()));
                                let seq_stats = utils::whole_seq_stats(&res);
                                let four_deg_codons =
                                    utils::gc_four_fold_deg_sites(&res, &degeneracy);
                                let gc_four = utils::four_fold_site_stats(four_deg_codons);

                                s.send(Output {
                                    fasta_id: fasta_record.id().to_string(),
                                    attr_id: attr_1.get("ID").unwrap().clone(),
                                    start: min_start,
                                    end: max_end,
                                    stats: seq_stats,
                                    gc_four_stats: gc_four,
                                    gc_3_stats,
                                })
                                .expect("Did not send.");
                                // clear the vector
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
                                let min_start: usize = *spliced
                                    .iter()
                                    .map(|x| x.start)
                                    .collect::<Vec<_>>()
                                    .first()
                                    .unwrap();
                                let max_end: usize = *spliced
                                    .iter()
                                    .map(|x| x.end)
                                    .collect::<Vec<_>>()
                                    .last()
                                    .unwrap();

                                let res = seq.join("");

                                let gc_3_stats =
                                    utils::gc_3(&res, Some(attr_1.get("ID").unwrap().clone()));
                                let seq_stats = utils::whole_seq_stats(&res);
                                let four_deg_codons =
                                    utils::gc_four_fold_deg_sites(&res, &degeneracy);
                                let gc_four = utils::four_fold_site_stats(four_deg_codons);

                                s.send(Output {
                                    fasta_id: fasta_record.id().to_string(),
                                    attr_id: attr_1.get("ID").unwrap().clone(),
                                    start: min_start,
                                    end: max_end,
                                    stats: seq_stats,
                                    gc_four_stats: gc_four,
                                    gc_3_stats,
                                })
                                .expect("Did not send.");

                                spliced.clear();
                            }
                        }
                    }
                }
            }
        });
    // collect and write.
    let res: Vec<Output> = receiver.iter().collect();
    let writer = Writer { output: res };

    eprintln!("[+]\tCollected output and writing.");

    writer.write(f)?;

    Ok(())
}
