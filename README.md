# GFF(3) stats

Given a genome and a corresponding GFF3 file, calculate various statistics on the coding regions.

Current output is a tsv, with bed-like first three columns (i.e. sequence ID, attribute ID, start, end...).

GC percent, GC skew, AT percent, and AC skew are calculated for each:
- raw CDS (or spliced CDS)
- four(six)-fold degenerate sites from the CDS
- third codon position for each codon in the CDS

gff-stats can also extract CDS/spliced CDS to a fasta file (see below).

The script takes into account strandedness and frame of the CDS. Performance is okay, but needs more testing on big genomes.

## Build

Building <a href="https://www.rust-lang.org/tools/install">requires Rust</a>. 

```bash
git clone https://github.com/tolkit/gff-stats
cd gff-stats
cargo build --release
# ./target/release/gff-stats is the executable
```

## Usage

`gff-stats -h`:

```bash
GFF(3) stats 0.2.0
Max Brown <mb39@sanger.ac.uk>
Extract GFF3 regions from a reference fasta and compute statistics on them.

USAGE:
    gff-stats [SUBCOMMAND]

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

SUBCOMMANDS:
    help    Prints this message or the help of the given subcommand(s)
    seq     Extract CDS regions to fasta format
    stat    Compute statistics on CDS regions
```

`gff-stats stat -h`

```bash
gff-stats-stat 
Compute statistics on CDS regions

USAGE:
    gff-stats stat [FLAGS] [OPTIONS] --fasta <fasta> --gff <gff>

FLAGS:
    -h, --help       Prints help information
    -p, --spliced    Compute stats on spliced CDS sequences?
    -V, --version    Prints version information

OPTIONS:
    -d, --degeneracy <degeneracy>    Calculate statistics on four-fold or six-fold (in addition to four-fold) degenerate
                                     codon sites. [default: fourfold]  [values: fourfold, sixfold]
    -f, --fasta <fasta>              The reference fasta file.
    -g, --gff <gff>                  The input gff file.
    -o, --output <output>            Output filename for the TSV (without extension). [default: gff-stat]
```

`gff-stats seq -h` 

```bash 
gff-stats-seq 
Extract CDS regions to fasta format

USAGE:
    gff-stats seq [FLAGS] [OPTIONS] --fasta <fasta> --gff <gff>

FLAGS:
    -h, --help       Prints help information
    -p, --spliced    Save the spliced extracted CDS fasta sequences?
    -V, --version    Prints version information

OPTIONS:
    -f, --fasta <fasta>      The reference fasta file.
    -g, --gff <gff>          The input gff file.
    -o, --output <output>    Output filename for the fasta (without extension). [default: gff-stat]
```

## TODO's

- Extract amino acid sequences