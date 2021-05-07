# GFF(3) stats

Given a genome and a corresponding GFF3 file, calculate various statistics on the coding regions.

Current output is a tsv, with bed-like first three columns (i.e. sequence ID, start, end...).

GC percent, GC skew, AT percent, and AC skew are calculated for each:
- raw CDS
- four(six)-fold degenerate sites from the CDS
- third codon position for each codon in the CDS

The script takes into account strandedness and frame of the CDS.

```bash
GFF stats 0.1.2
Max Brown <mb39@sanger.ac.uk>
Extract GFF3 regions from a reference fasta and compute statistics on them.

USAGE:
    gff-stats [OPTIONS] --fasta <fasta> --gff <gff>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -d, --degeneracy <degeneracy>    Calculate statistics on four-fold or six-fold (in addition to four-fold) degenerate
                                     codon sites. [default: fourfold]  [values: fourfold, sixfold]
    -f, --fasta <fasta>              The reference fasta file.
    -g, --gff <gff>                  The input gff file.
    -s, --sequences <sequences>      Save the extracted CDS fasta sequences? [default: false]
```