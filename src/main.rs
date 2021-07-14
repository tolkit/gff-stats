// gff-stats by Max Brown <mb39@sanger.ac.uk> //

use clap::{App, Arg, SubCommand};
use std::process;

use gff_stats::seq::seq;
use gff_stats::stat::stat;

// TODO: maybe return translated sequences?

fn main() {
    // command line options
    let matches = App::new("GFF(3) stats")
        .version(clap::crate_version!())
        .author("Max Brown <mb39@sanger.ac.uk>")
        .about("Extract GFF3 regions from a reference fasta and compute statistics on them.")
        .subcommand(
            SubCommand::with_name("stat")
            .about("Compute statistics on CDS regions")
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
            Arg::with_name("degeneracy")
                .short("d")
                .long("degeneracy")
                .required(false)
                .default_value("fourfold")
                .possible_values(&["fourfold", "sixfold"])
                .help("Calculate statistics on four-fold or six-fold (in addition to four-fold) degenerate codon sites."),
            )
        .arg(
            Arg::with_name("spliced")
                .short("p")
                .long("spliced")
                .help("Compute stats on spliced CDS sequences?"),
            )
            .arg(
                    Arg::with_name("output")
                        .short("o")
                        .long("output")
                        .takes_value(true)
                        .default_value("gff-stat")
                        .help("Output filename for the TSV (without extension)."),
                )
        )
        .subcommand(SubCommand::with_name("seq")
            .about("Extract CDS regions to fasta format")
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
                Arg::with_name("spliced")
                    .short("p")
                    .long("spliced")
                    .help("Save the spliced extracted CDS fasta sequences?"),
                )
                .arg(
                    Arg::with_name("output")
                        .short("o")
                        .long("output")
                        .takes_value(true)
                        .default_value("gff-stat")
                        .help("Output filename for the fasta (without extension)."),
                )
            )
            .get_matches();

    let subcommand = matches.subcommand();

    match subcommand.0 {
        "stat" => {
            let matches = subcommand.1.unwrap();
            stat::calculate_stats(matches);
        }
        "seq" => {
            let matches = subcommand.1.unwrap();
            seq::generate_seqs(matches);
        }
        _ => {
            eprintln!("Subcommand invalid, run with '--help' for subcommand options. Exiting.");
            process::exit(1);
        }
    }
}
