# GFF(3) stats

Given a genome and a corresponding GFF3 file, calculate various statistics on the coding regions (or extract them).

Current output is a tsv, with bed-like first three columns (i.e. sequence ID, attribute Parent ID, start, end...).

GC percent, GC skew, AT percent, and AC skew are calculated for each:
- raw CDS (or spliced CDS) (GC)
- four(six)-fold degenerate sites from the CDS (GC4)
- third codon position for each codon in the CDS (GC3)

gff-stats can also extract CDS/spliced CDS as a nucleotide or protein string to a fasta file (see below). Note this functionality is also provided by <a href="https://github.com/gpertea/gffread">gffread</a>.

Performance is okay, but needs more testing on big genomes.

Note that (obviously) the GFF must be in the correct format.

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
GFF(3) stats 0.2.1
Max Brown <mb39@sanger.ac.uk>
Extract GFF3 regions from a reference fasta and compute statistics on them.

USAGE:
    gff-stats [SUBCOMMAND]

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

SUBCOMMANDS:
    help    Prints this message or the help of the given subcommand(s)
    seq     Extract CDS regions to fasta format. Printed to stdout.
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
Extract CDS regions to fasta format. Printed to stdout.

USAGE:
    gff-stats seq [FLAGS] [OPTIONS] --fasta <fasta> --gff <gff>

FLAGS:
    -h, --help       Prints help information
    -p, --protein    Save the extracted CDS fasta sequences as a translated protein?
    -s, --spliced    Save the spliced extracted CDS fasta sequences?
    -V, --version    Prints version information

OPTIONS:
    -f, --fasta <fasta>      The reference fasta file.
    -g, --gff <gff>          The input gff file.
    -o, --output <output>    Output filename for the fasta (without extension). [default: gff-stat]
```