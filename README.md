# GFF(3) stats

Given a genome and a corresponding GFF3 file, calculate various statistics on the coding regions (or extract them).

Current output is a tsv, with bed-like first four columns (i.e. sequence ID, attribute Parent ID, start, end...).

GC percent, GC skew, AT percent, and AC skew are calculated for each:
- raw CDS (or spliced CDS) (GC)
- four(six)-fold degenerate sites from the CDS (GC4)
- third codon position for each codon in the CDS (GC3)

gff-stats can also extract CDS/spliced CDS as a nucleotide or protein string to a fasta file (see below). Note this functionality is also provided by <a href="https://github.com/gpertea/gffread">gffread</a> (see below). `gffread` may be faster as it indexes the fasta for quick random access.

## Build

Building <a href="https://www.rust-lang.org/tools/install">requires Rust</a>. 

```bash
git clone https://github.com/tolkit/gff-stats
cd gff-stats
cargo build --release
# ./target/release/gff-stats is the executable
# or
cargo install --path .
# to put gff-stats in your path
```

## Usage

### `gff-stats -h`

```bash
GFF(3) stats 0.2.2
Max Brown <mb39@sanger.ac.uk>
Extract GFF3 regions from a reference fasta and compute statistics on them.

USAGE:
    gff-stats [SUBCOMMAND]

OPTIONS:
    -h, --help       Print help information
    -V, --version    Print version information

SUBCOMMANDS:
    help    Print this message or the help of the given subcommand(s)
    seq     Extract CDS regions to fasta format. Printed to stdout.
    stat    Compute statistics on CDS regions
```

### `gff-stats stat -h`

```bash
gff-stats-stat 0.2.2
Compute statistics on CDS regions

USAGE:
    gff-stats stat [OPTIONS] --gff <gff> --fasta <fasta>

OPTIONS:
    -d, --degeneracy <degeneracy>    Calculate statistics on four-fold or six-fold (in addition to
                                     four-fold) degenerate codon sites. [default: fourfold]
                                     [possible values: fourfold, sixfold]
    -f, --fasta <fasta>              The reference fasta file.
    -g, --gff <gff>                  The input gff file.
    -h, --help                       Print help information
    -o, --output <output>            Output filename for the TSV (without extension). [default: gff-
                                     stat]
    -p, --spliced                    Compute stats on spliced CDS sequences?
    -V, --version                    Print version information
```

### `gff-stats seq -h` 

Cross testing with `gffread`:

```bash
# -x outputs spliced fastas
gffread -g ./tests/test_fasta.fna ./tests/test_gff.gff -x ./tests/test_gffread_x.fa
# equivalent to:
gff-stats seq -f ./tests/test_fasta.fna -g ./tests/test_gff.gff -s
# -y outputs spliced protein fastas
gffread -g ./tests/test_fasta.fna ./tests/test_gff.gff -y ./tests/test_gffread_y.fa
# equivalent to:
gff-stats seq -f ./tests/test_fasta.fna -g ./tests/test_gff.gff -sp
```

```bash 
gff-stats-seq 0.2.2
Extract CDS regions to fasta format. Printed to stdout.

USAGE:
    gff-stats seq [OPTIONS] --gff <gff> --fasta <fasta>

OPTIONS:
    -f, --fasta <fasta>      The reference fasta file.
    -g, --gff <gff>          The input gff file.
    -h, --help               Print help information
    -o, --output <output>    Output filename for the fasta (without extension). [default: gff-stat]
    -p, --protein            Save the extracted CDS fasta sequences as a translated protein?
    -s, --spliced            Save the spliced extracted CDS fasta sequences?
    -V, --version            Print version information
```

### Docs

<p align="center">
    <b>
        <a href="https://tolkit.github.io/gff-stats/">API documentation</a>
    </b>
</p>