use anyhow::Result;
use clap::{arg, crate_version, value_parser, ArgMatches, Command};
use std::path::PathBuf;

use gff_stats::seq;
use gff_stats::stat;
use gff_stats::Degeneracy;

fn main() -> Result<()> {
    // command line options
    let matches = Command::new("GFF(3) stats")
        .bin_name("gff-stats")
        .arg_required_else_help(true)
        .version(crate_version!())
        .author("Max Brown <mb39@sanger.ac.uk>")
        .about("Extract GFF3 regions from a reference fasta and compute statistics on them.")
        .subcommand(
            Command::new("stat")
                .about("Compute statistics on CDS regions.")
                .arg_required_else_help(true)
                .arg(
                    arg!(-g --gff <GFF> "The input GFF file")
                        .required(true)
                        .value_parser(value_parser!(PathBuf)),
                )
                .arg(arg!(-f --fasta <FASTA> "The reference fasta file")
                    .required(true)
                    .value_parser(value_parser!(PathBuf))
                )
                .arg(
                    arg!(-d --degeneracy [DEGENERACY] "Calculate statistics on four-fold or six-fold (in addition to four-fold) degenerate codon sites.")
                        .default_value("fourfold")
                        .value_parser(value_parser!(Degeneracy)),
                )
                .arg(
                    arg!(-s --spliced "Compute stats on spliced CDS sequences?")
                        .action(clap::ArgAction::SetTrue)
                )
                .arg(
                    arg!(-o --output [FILE])
                        .default_value("gff-stat")
                        .help("Output filename for the TSV (without extension)."),
                ),
        )
        .subcommand(
            Command::new("seq")
                .about("Extract CDS regions to fasta format. Printed to stdout.")
                .arg(
                    arg!(-g --gff <GFF> "The input GFF file")
                        .required(true)
                        .value_parser(value_parser!(PathBuf))
                )
                .arg(arg!(-f --fasta <FASTA> "The reference fasta file")
                    .required(true)
                    .value_parser(value_parser!(PathBuf))
                )
                .arg(
                    arg!(-s --spliced "Compute stats on spliced CDS sequences?")
                        .action(clap::ArgAction::SetTrue)
                )
                .arg(
                    arg!(-p --protein "Save the extracted CDS fasta sequences as a translated protein?")
                        .action(clap::ArgAction::SetTrue)
                )
                .arg(
                    arg!(-o --output [FILE])
                        .default_value("gff-stat")
                        .help("Output filename for the fasta (without extension)."),
                ),
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
        _ => unreachable!("Should never reach here."),
    }

    Ok(())
}
