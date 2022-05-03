use anyhow::Result;
use clap::{Arg, Command};
use std::process;

use gff_stats::seq;
use gff_stats::stat;

fn main() -> Result<()> {
    // command line options
    let matches = Command::new("GFF(3) stats")
        .version(clap::crate_version!())
        .propagate_version(true)
        .arg_required_else_help(true)
        .author("Max Brown <mb39@sanger.ac.uk>")
        .about("Extract GFF3 regions from a reference fasta and compute statistics on them.")
        .subcommand(
            Command::new("stat")
            .about("Compute statistics on CDS regions")
        .arg(
            Arg::new("gff")
            .short('g')
            .long("gff")
            .takes_value(true)
            .required(true)
            .help("The input gff file."),
        )
        .arg(
            Arg::new("fasta")
                .short('f')
                .long("fasta")
                .takes_value(true)
                .required(true)
                .help("The reference fasta file."),
        )
        .arg(
            Arg::new("degeneracy")
                .short('d')
                .long("degeneracy")
                .required(false)
                .default_value("fourfold")
                .possible_values(&["fourfold", "sixfold"])
                .help("Calculate statistics on four-fold or six-fold (in addition to four-fold) degenerate codon sites."),
            )
        .arg(
            Arg::new("spliced")
                .short('p')
                .long("spliced")
                .help("Compute stats on spliced CDS sequences?"),
            )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .takes_value(true)
                .default_value("gff-stat")
                .help("Output filename for the TSV (without extension)."),
            )
        )
        .subcommand(Command::new("seq")
            .about("Extract CDS regions to fasta format. Printed to stdout.")
            .arg(
            Arg::new("gff")
            .short('g')
            .long("gff")
            .takes_value(true)
            .required(true)
            .help("The input gff file."),
        )
        .arg(
            Arg::new("fasta")
                .short('f')
                .long("fasta")
                .takes_value(true)
                .required(true)
                .help("The reference fasta file."),
        )
        .arg(
            Arg::new("spliced")
                .short('s')
                .long("spliced")
                .help("Save the spliced extracted CDS fasta sequences?"),
        )
        .arg(
            Arg::new("protein")
                .short('p')
                .long("protein")
                .help("Save the extracted CDS fasta sequences as a translated protein?"),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .takes_value(true)
                .default_value("gff-stat")
                .help("Output filename for the fasta (without extension)."),
        )
    )
    .get_matches();

    let subcommand = matches.subcommand();

    match subcommand {
        Some(("stat", matches)) => {
            stat::calculate_stats(matches)?;
        }
        Some(("seq", matches)) => {
            seq::generate_seqs(matches)?;
        }
        _ => {
            eprintln!("Subcommand invalid, run with '--help' for subcommand options. Exiting.");
            process::exit(1);
        }
    }

    Ok(())
}
